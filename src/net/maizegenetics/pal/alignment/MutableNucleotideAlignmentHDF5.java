/*
 * MutableBitNucleotideAlignmentHDF5
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import org.apache.log4j.Logger;

/**
 *
 * @author terry & ed
 */
public class MutableNucleotideAlignmentHDF5 extends AbstractAlignment implements MutableAlignment {

    private static final Logger myLogger = Logger.getLogger(MutableNucleotideAlignmentHDF5.class);
    private boolean myIsDirty = true;
    private byte[][] myData;
    private List<Identifier> myIdentifiers = new ArrayList<Identifier>();
//    private final int myMaxTaxa;
    private final int myMaxNumSites;
    private int myNumSites = 0;
    private int myNumSitesStagedToRemove = 0;
    protected int[] myVariableSites;
    private List<Locus> myLocusToLociIndex = new ArrayList<Locus>();
//    protected int[] myLocusIndices;
    private int[] myLocusOffsets = null;
    protected String[] mySNPIDs;
    IHDF5Writer myWriter=null;
    private HashMap<Integer,byte[]> myDataHashMap;
    private byte[] singleCache=null;
    private int taxonCached=-1;
    HDF5IntStorageFeatures genoFeatures = HDF5IntStorageFeatures.createDeflation(HDF5IntStorageFeatures.MAX_DEFLATION_LEVEL);

    protected MutableNucleotideAlignmentHDF5(Alignment a, int maxNumTaxa, int maxNumSites) {
        super(a.getAlleleEncodings());
        myMaxNumAlleles = a.getMaxNumAlleles();

        if (a.getAlleleEncodings().length != 1) {
            throw new IllegalArgumentException("MutableBitNucleotideAlignmentHDF5: init: must only have one allele encoding.");
        }

        if (a.getSiteCount() > maxNumSites) {
            throw new IllegalArgumentException("MutableBitNucleotideAlignmentHDF5: init: initial number of sites can't be more than max number of sites.");
        }

        if (a.getSequenceCount() > maxNumTaxa) {
            throw new IllegalArgumentException("MutableBitNucleotideAlignmentHDF5: init: initial number of taxa can't be more than max number of taxa.");
        }
        myMaxNumSites = maxNumSites;
        myNumSites = a.getSiteCount();
        myReference = new byte[myMaxNumSites];
        Arrays.fill(myReference, Alignment.UNKNOWN_DIPLOID_ALLELE);
        initData();
        initTaxa(a.getIdGroup());
        loadAlleles(a);
        //loadLoci(a);
        System.arraycopy(a.getSNPIDs(), 0, mySNPIDs, 0, a.getSiteCount());
        System.arraycopy(a.getPhysicalPositions(), 0, myVariableSites, 0, a.getSiteCount());
    }

    public static MutableNucleotideAlignmentHDF5 getInstance(Alignment a, int maxTaxa, int maxNumSites) {
        return MutableNucleotideAlignmentHDF5.getInstance(a, maxTaxa, maxNumSites);
    }

    public static MutableNucleotideAlignmentHDF5 getInstance(Alignment a) {
        return getInstance(a, a.getSequenceCount(), a.getSiteCount());
    }


    protected MutableNucleotideAlignmentHDF5(IHDF5Writer reader, List<Identifier> idGroup, int[] variableSites, 
            List<Locus> locusToLociIndex, int[] lociOffsets, String[] siteNames) {
        super(NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);

        if (variableSites.length != siteNames.length) {
            throw new IllegalArgumentException("MutableBitNucleotideAlignmentHDF5: init: number variable sites, loci, and site names must be same.");
        }
        myWriter=reader;
        myMaxNumSites = siteNames.length;
        myNumSites = siteNames.length;

        myVariableSites = variableSites;
        myLocusToLociIndex = locusToLociIndex;
        myLocusOffsets=lociOffsets;
      //  myLocusIndices = locusIndices;
        mySNPIDs = siteNames;

        myReference = new byte[myMaxNumSites];
        Arrays.fill(myReference, Alignment.UNKNOWN_DIPLOID_ALLELE);
        Collections.sort(idGroup);
       // idGroup
        
        myIdentifiers = idGroup;
    }

//    public static MutableNucleotideAlignmentHDF5 getInstance(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
//        return new MutableNucleotideAlignmentHDF5(idGroup, initNumSites, maxNumTaxa, maxNumSites);
//    }
//
//    public static MutableNucleotideAlignmentHDF5 getInstance(IdGroup idGroup, int initNumSites) {
//        return MutableNucleotideAlignmentHDF5.getInstance(idGroup, initNumSites, 
//                idGroup.getIdCount(), initNumSites);
//    }

    public static MutableNucleotideAlignmentHDF5 getInstance(String filename) {
        IHDF5Writer reader = HDF5Factory.open(filename);
        
        //derive the taxa list on the fly 
        List<HDF5LinkInformation> fields=reader.getAllGroupMemberInformation(HapMapHDF5Constants.GENOTYPES, true);
        ArrayList<Identifier> taxaList=new ArrayList<Identifier>();
        for (HDF5LinkInformation is : fields) {
            if(is.isDataSet()==false) continue;
   //         System.out.println(is.getName());
            taxaList.add(new Identifier(is.getName()));
        }

 //       IdGroup idgroup = new SimpleIdGroup(taxa);

        byte[][] alleles = reader.readByteMatrix(HapMapHDF5Constants.ALLELES);

        MDArray<String> alleleStatesMDArray = reader.readStringMDArray(HapMapHDF5Constants.ALLELE_STATES);
        int[] dimensions = alleleStatesMDArray.dimensions();
        int numEncodings = dimensions[0];
        int numStates = dimensions[1];
        String[][] alleleStates = new String[numEncodings][numStates];
        for (int e = 0; e < numEncodings; e++) {
            for (int s = 0; s < numStates; s++) {
                alleleStates[e][s] = alleleStatesMDArray.get(e, s);
            }
        }

        int[] variableSites = reader.readIntArray(HapMapHDF5Constants.POSITIONS);
        int maxNumAlleles = reader.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.MAX_NUM_ALLELES);
        boolean retainRare = reader.getBooleanAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.RETAIN_RARE_ALLELES);
        String[] lociStrings = reader.readStringArray(HapMapHDF5Constants.LOCI);
        int numLoci = lociStrings.length;
       // Locus[] loci = new Locus[numLoci];
        ArrayList<Locus> loci=new ArrayList<Locus>();
        for (String lS : lociStrings) {loci.add(new Locus(lS));}
        int[] lociOffsets = reader.readIntArray(HapMapHDF5Constants.LOCUS_OFFSETS);
        String[] snpIds = reader.readStringArray(HapMapHDF5Constants.SNP_IDS);
        return new MutableNucleotideAlignmentHDF5(reader, taxaList, variableSites, loci, lociOffsets, snpIds);
    }
    
    private void initData() {
        myVariableSites = new int[myMaxNumSites];
        Arrays.fill(myVariableSites, -1);
        mySNPIDs = new String[myMaxNumSites];
        Arrays.fill(mySNPIDs, null);
    }

    private void initTaxa(IdGroup idGroup) {
        for (int i = 0, n = idGroup.getIdCount(); i < n; i++) {
            myIdentifiers.add(idGroup.getIdentifier(i));
        }
    }

    private void loadAlleles(Alignment a) {
        int numSites = a.getSiteCount();
        int numSeqs = a.getSequenceCount();

        for (int s = 0; s < numSites; s++) {
            for (int t = 0; t < numSeqs; t++) {
                myData[t][s] = a.getBase(t, s);
            }
        }
    }

    public byte getBase(int taxon, int site) {
        if(taxon!=taxonCached) {
            singleCache=myWriter.readByteArray(HapMapHDF5Constants.GENOTYPES + "/" + getFullTaxaName(taxon));
            taxonCached=taxon;
        }
        return singleCache[site];
       // return myData[taxon][site];
        
        //http://stackoverflow.com/questions/12319741/limited-size-hash-map
    }

    public boolean isSBitFriendly() {
        return false;
    }

    public boolean isTBitFriendly() {
        return false;
    }

    @Override
    public int getSiteCount() {
        return myNumSites;
    }

    @Override
    public int getSequenceCount() {
        return myIdentifiers.size();
    }

    @Override
    public IdGroup getIdGroup() {
        Identifier[] ids = new Identifier[myIdentifiers.size()];
        myIdentifiers.toArray(ids);
        return new SimpleIdGroup(ids);
    }

    @Override
    public String getTaxaName(int index) {
        return myIdentifiers.get(index).getName();
    }

    @Override
    public String getFullTaxaName(int index) {
        return myIdentifiers.get(index).getFullName();
    }

    @Override
    public boolean isPhased() {
        return false;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return true;
    }

    @Override
    public String getGenomeAssembly() {
        return "AGPV2";
    }

    @Override
    public int[] getPhysicalPositions() {
        return myVariableSites.clone();
    }

    @Override
    public int getPositionInLocus(int site) {
        try {
            if (myVariableSites[site] < 0) {
                return site;
            }
            return myVariableSites[site];
        } catch (Exception e) {
            return site;
        }
    }

//    @Override
//    public Locus getLocus(int site) {
//        return myLocusToLociIndex.get(myLocusIndices[site]);
//    }

    @Override
    public Locus[] getLoci() {
        Locus[] result = new Locus[myLocusToLociIndex.size()];
        for (int i = 0; i < myLocusToLociIndex.size(); i++) {
            result[i] = myLocusToLociIndex.get(i);
        }
        return result;
    }

    @Override
    public int getNumLoci() {
        return myLocusToLociIndex.size();
    }



    @Override
    public int[] getStartAndEndOfLocus(Locus locus) {

        if (isDirty()) {
            throw new IllegalStateException("MutableBitNucleotideAlignmentHDF5: getStartAndEndOfLocus: this alignment is dirty.");
        }

        Locus[] loci = getLoci();
        int[] lociOffsets = getLociOffsets();
        int numLoci = getNumLoci();
        for (int i = 0; i < numLoci; i++) {
            if (locus.equals(loci[i])) {
                int end = 0;
                if (i == numLoci - 1) {
                    end = getSiteCount();
                } else {
                    end = lociOffsets[i + 1];
                }
                return new int[]{lociOffsets[i], end};
            }
        }
        throw new IllegalArgumentException("AbstractAlignment: getStartAndEndOfLocus: this locus not defined: " + locus.getName());
    }

    @Override
    public String[] getSNPIDs() {
        return mySNPIDs;
    }

    @Override
    public String getSNPID(int site) {
        if ((mySNPIDs == null) || (mySNPIDs.length == 0) || (mySNPIDs[site] == null)) {
            return "S" + getLocus(site).getChromosomeName() + "_" + getPositionInLocus(site);
        }
        return mySNPIDs[site];
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {

        if (isDirty()) {
            throw new IllegalStateException("MutableBitNucleotideAlignmentHDF5: getSiteOfPhysicalPosition: this alignment is dirty.");
        }

        if (myVariableSites == null) {
            return physicalPosition;
        }
        try {
            if (locus == null) {
                locus = myLocusToLociIndex.get(0);
            }
            int[] startEnd = getStartAndEndOfLocus(locus);
            return Arrays.binarySearch(myVariableSites, startEnd[0], startEnd[1], physicalPosition);
        } catch (Exception e) {
            e.printStackTrace();
            return -1;
        }
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus, String snpID) {

        if (isDirty()) {
            throw new IllegalStateException("MutableBitNucleotideAlignmentHDF5: getSiteOfPhysicalPosition: this alignment is dirty.");
        }

        if (myVariableSites == null) {
            return physicalPosition;
        }
        if (locus == null) {
            locus = myLocusToLociIndex.get(0);
        }
        int[] startEnd = getStartAndEndOfLocus(locus);
        int result = Arrays.binarySearch(myVariableSites, startEnd[0], startEnd[1], physicalPosition);
        if (result < 0) {
            return result;
        } else {
            if (snpID.equals(getSNPID(result))) {
                return result;
            } else {
                int index = result - 1;
                while ((index >= startEnd[0]) && (getPositionInLocus(index) == physicalPosition)) {
                    if (snpID.equals(getSNPID(index))) {
                        return index;
                    }
                    index--;
                }
                index = result + 1;
                while ((index < startEnd[1]) && (getPositionInLocus(index) == physicalPosition)) {
                    if (snpID.equals(getSNPID(index))) {
                        return index;
                    }
                    index++;
                }
                return -result - 1;
            }
        }

    }

    @Override
    public boolean retainsRareAlleles() {
        return false;
    }

    @Override
    public byte[] getAlleles(int site) {
        //TODO this should be just stored in the HDF5 file
        //reprocess as needed if taxa added
        int[][] alleles = getAllelesSortedByFrequency(site);
        int resultSize = alleles[0].length;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i];
        }
        return result;
    }

    // Mutable Methods...
    public void setBase(int taxon, int site, byte newBase) {
        myData[taxon][site] = newBase;
    }

    public void setBase(Identifier taxon, String siteName, Locus locus, int physicalPosition, byte newBase) {

        int taxonIndex = myIdentifiers.indexOf(taxon);
        if (taxonIndex == -1) {
            throw new IllegalArgumentException("MutableBitNucleotideAlignmentHDF5: setBase: taxon not found.");
        }

        int site = getSiteOfPhysicalPosition(physicalPosition, locus, siteName);
        if (site < 0) {
            throw new IllegalStateException("MutableBitNucleotideAlignmentHDF5: setBase: physical position: " + physicalPosition + " in locus: " + locus.getName() + " not found.");
        } else {
            if (!siteName.equals(getSNPID(site))) {
                throw new IllegalStateException("MutableBitNucleotideAlignmentHDF5: setBase: site names at physical position: " + physicalPosition + " in locus: " + locus.getName() + " does not match: " + siteName);
            }
        }

        myData[taxonIndex][site] = newBase;

    }

    public void setBaseRange(int taxon, int startSite, byte[] newBases) {
        for (int i = 0; i < newBases.length; i++) {
            myData[taxon][startSite++] = newBases[i];
        }
    }

    public void setReferenceAllele(int site, byte diploidAllele) {
        myReference[site] = diploidAllele;
    }

    public void addTaxon(Identifier id) {
        String basesPath = HapMapHDF5Constants.GENOTYPES + "/" + id.getFullName();
        myWriter.createByteArray(basesPath, myNumSites, genoFeatures);
        myIdentifiers.add(id);
    }

    public void setTaxonName(int taxon, Identifier id) {
        if (taxon >= myIdentifiers.size()) {
            throw new IllegalStateException("MutableBitNucleotideAlignmentHDF5: setTaxonName: this taxa index does not exist: " + taxon);
        }
        String currentPath = HapMapHDF5Constants.GENOTYPES + "/" + myIdentifiers.get(taxon);
        String newPath = HapMapHDF5Constants.GENOTYPES + "/" + id.getFullName();
        myWriter.move(currentPath, newPath);
        myIdentifiers.set(taxon, id);
    }

    public void removeTaxon(int taxon) {
        myIdentifiers.remove(taxon);
        String currentPath = HapMapHDF5Constants.GENOTYPES + "/" + myIdentifiers.get(taxon);
        myWriter.delete(currentPath);
        //TODO remove from cache

    }

    public void clean() {
        myIsDirty = false;
        myNumSites -= myNumSitesStagedToRemove;
        myNumSitesStagedToRemove = 0;
    }

    public boolean isDirty() {
        return myIsDirty;
    }

    private void setDirty() {
        myLocusOffsets = null;
        myIsDirty = true;
    }

    private void setClean() {
        myIsDirty = false;
    }



    private int getLocusIndex(Locus locus) {
        for (int i = 0; i < myLocusToLociIndex.size(); i++) {
            if (myLocusToLociIndex.get(i).equals(locus)) {
                return i;
            }
        }
        return -1;
    }

    public void setDepthForAlleles(int taxon, int site, byte[] values) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public void setCommonAlleles(int site, byte[] values) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void addSite(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void removeSite(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void clearSiteForRemoval(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setPositionOfSite(int site, int position) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setLocusOfSite(int site, Locus locus) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
