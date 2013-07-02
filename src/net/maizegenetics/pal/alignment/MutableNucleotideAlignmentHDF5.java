/*
 * MutableBitNucleotideAlignmentHDF5
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import static net.maizegenetics.pal.alignment.Alignment.UNKNOWN_DIPLOID_ALLELE;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.ProgressListener;
import org.apache.log4j.Logger;


/**
 * MutableNucleotideAlignmentHDF5 is a byte based HDF5 mutable alignment.  The mutability 
 * only allows for the addition and deletion of taxa, and for the modification of genotypes.
 * It does NOT permit the insertion or deletion of sites (although it does allow merging of
 * alignments).  Deletion may be possible in the future.
 * <p>
 * The class uses a local cache to provide performance.  The default cache is defined by 
 * defaultCacheSize, and is set 4096 chunks.  The genotypes in the cache are chunked based 
 * on defaultSiteCache also set to default of 4096 sites (4096 * 4096 = 16M sites).
 * The cache be filled with any balance of chunks, and the first added are the first 
 * removed if the cache fills.  If there will be lots of reading from the alignment, set the 
 * defaultCacheSize > the number of taxa.  If you will only be writing to the object,
 * then set the defaultCacheSize small.
 * <p>
 * This class is very speed efficient at most operations, but it is slow at setBase() because it 
 * does not cache the bases being written for thread safety.  Write genotypes in bulk with
 * setBaseRange().
 * <p>
 * To create this file format use:
 * ExportUtils.writeToMutableHDF5(Alignment a, String newHDF5file, boolean includeGenotypes)
 * To instantiate it:
 * MutableNucleotideAlignmentHDF5.
 * <p>
 * If you use this in a research publication - cite Bradbury et al
 * 
 * <p>
 * TODO:
 * <li>Add support for bit encoding - but just two states (Major/Minor) or (Ref/Alt)
 * <li>Add efficient support for merging chromosomes
 * 
 * @author Ed Buckler & Terry Casstevens
 */
public class MutableNucleotideAlignmentHDF5 extends AbstractAlignment implements MutableAlignment {

    private static final Logger myLogger = Logger.getLogger(MutableNucleotideAlignmentHDF5.class);
    private String fileName;
    private boolean myIsDirty = true;
    private List<Identifier> myIdentifiers = new ArrayList<Identifier>();
  
    
    private final int myMaxNumSites;  //try to delete
    private int myNumSitesStagedToRemove = 0;  //try to delete
    
    
    protected int[] myVariableSites;
    private List<Locus> myLocusList = new ArrayList<Locus>();
    protected int[] myLocusIndices;  //consider changing to byte or short to save memory, generally on 10 chromosomes
    private int[] myLociOffsets = null;

    //Site descriptor in memory
    private byte[][] myAlleleFreqOrder=null;
    private int[][] myAlleleCnt=null;
    private float[] maf=null;
    private float[] siteCov=null;
    protected String[] mySNPIDs;  //perhaps don't keep in memory, consider compressed
    
    IHDF5Writer myWriter=null;
    

    HDF5IntStorageFeatures genoFeatures = HDF5IntStorageFeatures.createDeflation(2);
    private int defaultCacheSize=1<<16;
    private int defaultSiteCache=1<<16;
    
    /*Byte representations of DNA sequences are stored in blocks of 65536 sites*/
    public static final int hdf5GenoBlock=1<<16;
    public static final int hdf5GenoBlockMask=~(hdf5GenoBlock-1);
//    public boolean favorCachingOfEntireTaxon=true;

    private LoadingCache<Long,byte[]> myGenoCache; //key (taxa <<< 32)+BLOCK
    private CacheLoader<Long,byte[]> genoLoader = new CacheLoader<Long,byte[]>() {
        public byte[] load(Long key) {
            int[] taxonSite=getTaxonSiteFromKey(key);
            long offset=taxonSite[1]<<16;
            byte[] data=new byte[0];
            synchronized(myWriter) {
                data=myWriter.readAsByteArrayBlockWithOffset(getTaxaGenoPath(taxonSite[0]),1<<16,offset);}
            return data;
    } };
        
    private boolean cacheDepth=false;
    private LRUCache<Long,byte[][]> myDepthCache=null;  //key (taxa <<< 32)+startSite
    
  
    private LoadingCache<Integer,BitSet[]> tBitCache; //taxon id number, Bitset[2] = {MajorBitSet, MinorBitSet}    
    private CacheLoader<Integer,BitSet[]> tBitLoader = new CacheLoader<Integer,BitSet[]>() {
        public BitSet[] load(Integer key) {
            BitSet[] bs=AlignmentUtils.calcBitPresenceFromGenotype(getBaseRow(key), myAlleleFreqOrder[0], myAlleleFreqOrder[1]);
            return bs;
    } };
    
    private IdGroup myIdGroup;  //set to null whenever dirty

    protected MutableNucleotideAlignmentHDF5(String fileName, IHDF5Writer reader, List<Identifier> idGroup, int[] variableSites, 
            List<Locus> locusList, int[] locusIndices, String[] siteNames, int defaultCacheSize) {
        super(NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        myMaxNumAlleles=NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES;
        this.fileName=fileName;
        this.defaultCacheSize=defaultCacheSize;
        
        if (variableSites.length != siteNames.length) {
            throw new IllegalArgumentException("MutableBitNucleotideAlignmentHDF5: init: number variable sites, loci, and site names must be same.");
        }
        if (variableSites.length != siteNames.length) {
            throw new IllegalArgumentException("MutableBitNucleotideAlignmentHDF5: init: number variable sites, loci, and site names must be same.");
        }
        myWriter=reader;
        myMaxNumSites = siteNames.length;
        myNumSites = siteNames.length;
        this.defaultSiteCache=(myNumSites<defaultSiteCache)?myNumSites:defaultSiteCache;

        myVariableSites = variableSites;
        myLocusList = locusList;
        myLocusIndices = locusIndices;
        mySNPIDs = siteNames;

        myReference = new byte[myMaxNumSites];
        Arrays.fill(myReference, Alignment.UNKNOWN_DIPLOID_ALLELE);
        Collections.sort(idGroup);
        
        myIdentifiers = idGroup;
        recalcLocusCoordinates();
        initAllDataSupport();
        loadSiteDescriptorsToMemory();
        initGenotypeCache();
        myIsDirty=false;
    }

    public static MutableNucleotideAlignmentHDF5 getInstance(String filename) {
        return MutableNucleotideAlignmentHDF5.getInstance(filename, 1<<16);
    }
    
    public static boolean isMutableNucleotideAlignmentHDF5(String filename) {
        IHDF5Reader reader = HDF5Factory.open(filename);
        boolean result=reader.exists(HapMapHDF5Constants.GENOTYPES);
        reader.close();
        return result;
    }
    
    public static MutableNucleotideAlignmentHDF5 getInstance(String filename, int defaultCacheSize) {
        IHDF5Writer reader = HDF5Factory.open(filename);
        //derive the taxa list on the fly 
        List<HDF5LinkInformation> fields=reader.getAllGroupMemberInformation(HapMapHDF5Constants.GENOTYPES, true);
        ArrayList<Identifier> taxaList=new ArrayList<Identifier>();
        for (HDF5LinkInformation is : fields) {
            if(is.isDataSet()==false) continue;
            taxaList.add(new Identifier(is.getName()));
        }
        
        //TODO this is a good idea, it should be passed into the constructor and stored
//        byte[][] alleles = reader.readByteMatrix(HapMapHDF5Constants.ALLELES);

        int[] variableSites = reader.readIntArray(HapMapHDF5Constants.POSITIONS);
        String[] lociStrings = reader.readStringArray(HapMapHDF5Constants.LOCI);
        ArrayList<Locus> loci=new ArrayList<Locus>();
        for (String lS : lociStrings) {loci.add(new Locus(lS));}
        int[] locusIndices;
        if(reader.exists(HapMapHDF5Constants.LOCUS_INDICES)) {
            locusIndices = reader.readIntArray(HapMapHDF5Constants.LOCUS_INDICES);}
        else {locusIndices = new int[variableSites.length];  //this to support an old format.
            Arrays.fill(locusIndices, 0); //fill with zero
        }
        String[] snpIds = reader.readStringArray(HapMapHDF5Constants.SNP_IDS);
        return new MutableNucleotideAlignmentHDF5(filename, reader, taxaList, variableSites, loci, locusIndices, snpIds, defaultCacheSize);
    }
    
    
    private void initAllDataSupport() {
        if(!myWriter.exists(HapMapHDF5Constants.SBIT)) myWriter.createGroup(HapMapHDF5Constants.SBIT);
        if(!myWriter.exists(HapMapHDF5Constants.TBIT)) myWriter.createGroup(HapMapHDF5Constants.TBIT);
        if(!myWriter.exists(HapMapHDF5Constants.DEPTH)) myWriter.createGroup(HapMapHDF5Constants.DEPTH);
        if(!myWriter.exists(HapMapHDF5Constants.SITE_DESC)) clean();
    }

    
    private void initGenotypeCache() {
        int minLoad=(3*getSequenceCount())/2;
        myGenoCache= CacheBuilder.newBuilder()
            .maximumSize(minLoad)
            .build(genoLoader);
        tBitCache = CacheBuilder.newBuilder()
            .maximumSize(1024)
            .build(tBitLoader);
    }
    
    private void initDepthCache() {
        myDepthCache=new LRUCache<Long, byte[][]>(defaultCacheSize);
//        int minLoad=(defaultCacheSize<getSequenceCount())?defaultCacheSize:getSequenceCount();
//        int start=(0/defaultSiteCache)*defaultSiteCache;
//        int realSiteCache=(myNumSites-start<defaultSiteCache)?myNumSites-start:defaultSiteCache;
//     //   int xSizeCache=3000;
//        for (int i = 0; i<minLoad; i++) {
////            byte[][] test2=myWriter.readByteMatrix(getTaxaDepthPath(i));
////            System.out.println(Arrays.deepToString(test2));
////            byte[][] test=myWriter.readByteMatrixBlockWithOffset(getTaxaDepthPath(i), 6, xSizeCache, 0, start);
////            System.out.println(Arrays.deepToString(test));
//            myDepthCache.put(getCacheKey(i,0), myWriter.readByteMatrixBlockWithOffset(getTaxaDepthPath(i), 6, realSiteCache, 0, start));
//        }
    }
    
    protected static long getCacheKey(int taxon, int site) {
        return ((long)taxon<<32)+(site/hdf5GenoBlock);
    }
    
    protected static int[] getTaxonSiteFromKey(long key) {
        return new int[]{(int)(key>>>32),(int)((key<<32)>>>32)};
    }
    
    public void clearGenotypeCache() {
        //myGenoCache.clear();
        myGenoCache.invalidateAll();
    }
    


    
//    private byte[] cacheTaxonSiteBlock(int taxon, int site) {
//        return cacheTaxonSiteBlock(taxon, site, getCacheKey(taxon, site));
//    }
    
    private byte[][] cacheDepthBlock(int taxon, int site, long key) {
        int start=(site/defaultSiteCache)*defaultSiteCache;
        int realSiteCache=(myNumSites-start<defaultSiteCache)?myNumSites-start:defaultSiteCache;
        byte[][] data=myWriter.readByteMatrixBlockWithOffset(getTaxaDepthPath(taxon), 6, realSiteCache, 0, start);
        if(data==null) return null;
        myDepthCache.put(key, data);
        return data;
    }
    
    private byte[][] cacheDepthBlock(int taxon, int site) {
        return cacheDepthBlock(taxon, site, getCacheKey(taxon, site));
    }
    
    public int getCacheSize() {
        return (int)myGenoCache.size();
    }
    
//    private void removeCacheTaxonSiteBlock(int taxon, int site) {
//        long key=getCacheKey(taxon, site);
//        myDataCache.remove(key);
//    }
    
    protected String getTaxaGenoPath(int taxon) {
        return HapMapHDF5Constants.GENOTYPES + "/" + getFullTaxaName(taxon);
    }
    
    protected String getTaxaDepthPath(int taxon) {
        return HapMapHDF5Constants.DEPTH + "/" + getFullTaxaName(taxon);
    }

    @Override
    public byte getBase(int taxon, int site) {
        long key=getCacheKey(taxon,site);
        try{
            byte[] data=myGenoCache.get(key);
            return data[site%hdf5GenoBlock];
        } catch(ExecutionException e) {
           e.printStackTrace();
           throw new UnsupportedOperationException("Error in getBase from cache");
        }
        //if((taxon==100)&&(site==100)) System.out.println("Called Get in getBase");
//        if(data==null) {data=cacheTaxonSiteBlock(taxon, site, key);}
//        return data[site%hdf5GenoBlock];
    }
    
    @Override
    public byte[] getBaseRow(int taxon) {
        return getBaseRange(taxon, 0, myNumSites);
    }
    
    @Override
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {
        byte[] result = new byte[endSite - startSite];
        ByteBuffer bbR=ByteBuffer.wrap(result);
        for (int site = startSite; site < endSite; site+=hdf5GenoBlock) {
            long key=getCacheKey(taxon,site);
            byte[] data;
            try{data=myGenoCache.get(key);
            } catch(ExecutionException e) {
               e.printStackTrace();
               throw new UnsupportedOperationException("Error in getBaseRange from cache");
            }
            if(data==null) return null;
            int blockStart=site&hdf5GenoBlockMask;
            int offset=(blockStart<startSite)?(startSite-blockStart):0;
            int length=(endSite-site>hdf5GenoBlock)?data.length:(endSite-site);
            bbR.put(data, offset, length);
        }
        return result;
    }
    
    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon,site));
    }
    
   @Override
    public String getBaseAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            builder.append(getBaseAsString(taxon, i));
        }
        return builder.toString();
    }
    
    @Override
    public byte[] getDepthForAlleles(int taxon, int site) {
        if(myDepthCache==null) initDepthCache();
        long key=getCacheKey(taxon,site);
        byte[][] data=myDepthCache.get(key);
        if(data==null) {data=cacheDepthBlock(taxon, site, key);}
        byte[] result=new byte[6];
        for (int i = 0; i < 6; i++) {
            result[i]=data[i][site%defaultSiteCache];
        }
        return result;
    }
    
    

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        if(alleleNumber>1) throw new UnsupportedOperationException("Will not be supported.  Only top two allele reported");
        try {
            BitSet[] bs=tBitCache.get(taxon);
            return bs[alleleNumber];
        } catch (ExecutionException e) {
            e.printStackTrace();
        }
        return null;
    }
    

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        BitSet result=getAllelePresenceForAllSites(taxon, alleleNumber);
        return result.getBits(startBlock, endBlock-1);
    }

    @Override
    public void optimizeForTaxa(ProgressListener listener) {
        //do nothing now, later setup 
        
    }
    
    

    @Override
    public boolean isSBitFriendly() {
        return false;
    }

    @Override
    public boolean isTBitFriendly() {
        return true;
    }

    @Override
    public int getTotalNumAlleles() {
        return 6; 
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
        if(myIdGroup==null) {
            Identifier[] ids = new Identifier[myIdentifiers.size()];
            myIdentifiers.toArray(ids);
            myIdGroup=new SimpleIdGroup(ids);
        } 
        return myIdGroup;
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
    
    private void recalcLocusCoordinates() {
        myLociOffsets=new int[myLocusList.size()];
        int currStart=0;
        int currLIndex=myLocusIndices[currStart];
        myLociOffsets[0]=currStart;
        for (int s = 0; s < myNumSites; s++) {
            if(currLIndex!=myLocusIndices[s]) {
                String name=getLocus(currStart).getName();
                myLocusList.set(currLIndex, new Locus(name, name, currStart, s-1, null, null));
                currStart=s;
                currLIndex=myLocusIndices[s];
                myLociOffsets[currLIndex]=currStart;
            }   
        }
        String name=getLocus(currStart).getName();
        myLocusList.set(currLIndex, new Locus(name, name, currStart, myNumSites-1, null, null));
    }


    @Override
    public Locus[] getLoci() {
        Locus[] result = new Locus[myLocusList.size()];
        for (int i = 0; i < myLocusList.size(); i++) result[i] = myLocusList.get(i);
        return result;
    }

    @Override
    public int getNumLoci() {
        return myLocusList.size();
    }

    @Override
    public Locus getLocus(int site) {
        return myLocusList.get(myLocusIndices[site]);
    }
    
    @Override
    public int[] getLociOffsets() {
        return myLociOffsets;
    }
    

    @Override
    public int[] getStartAndEndOfLocus(Locus locus) {
        for (Locus mL : myLocusList) {
            if(mL.getName().equals(locus.getName())) return new int[]{mL.getStart(),mL.getEnd()};
        }
        return null;
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
                locus = myLocusList.get(0);
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
            locus = myLocusList.get(0);
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
        byte[] result=new byte[6];
        int i;
        for (i = 0; (myAlleleFreqOrder[i][site]!=Alignment.UNKNOWN_ALLELE)&&(i<6); i++) {
            result[i]=myAlleleFreqOrder[i][site];
        }
        return Arrays.copyOf(result, i);
    }
    
   @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        byte[] alleles=getAlleles(site);
        int result[][] = new int[2][alleles.length];
//        System.out.printf("gASBF AL:%d s:%d %d %d %d %d %n", alleles.length, site, myAlleleCnt[0][site], myAlleleCnt[1][site],  myAlleleCnt[2][site], myAlleleCnt[3][site]);
        for (int i = 0; i < alleles.length; i++) {
           result[0][i]=alleles[i];
            
           result[1][i]=myAlleleCnt[alleles[i]][site];
       }
       return result;
    }
    
    public void setCacheHints(int numberOfTaxa, int numberOfSites) {
        //Number of taxa to keep in the cache
        //size of the block of sequence to cache.
        //linkedHashList <Taxon, <startSite,byte[]>>
    }

    // Mutable Methods...
    public void setBase(int taxon, int site, byte newBase) {
        this.myWriter.writeByteArrayBlock(getTaxaGenoPath(taxon),new byte[]{newBase},site);
        //cacheTaxonSiteBlock(taxon, site);
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
        setBase(taxonIndex,site,newBase);
        myIsDirty=true;
        //myData[taxonIndex][site] = newBase;

    }

    public synchronized void setBaseRange(int taxon, int startSite, byte[] newBases) {
        myWriter.writeByteArrayBlock(getTaxaGenoPath(taxon),newBases,startSite);
        for (int site = startSite; site < (startSite+newBases.length); site+=defaultSiteCache) {
            //cacheTaxonSiteBlock(taxon, site);
        }
        myIsDirty=true;
        
    }
    
    public synchronized void setAllBases(int taxon, byte[] newBases) {
       myWriter.writeByteArray(getTaxaGenoPath(taxon), newBases, genoFeatures);
       for (int site = 0; site < newBases.length; site+=defaultSiteCache) {
//            cacheTaxonSiteBlock(taxon, site);
        }
       myIsDirty=true;
    }
    
    public synchronized void setAllDepth(int taxon, byte[][] depth) {
        if(depth.length!=6) throw new IllegalStateException("Just set A, C, G, T, -, + all at once");
        if(depth[0].length!=myNumSites) throw new IllegalStateException("Setting all depth.  Wrong number of sites");
        myWriter.writeByteMatrix(getTaxaDepthPath(taxon), depth, genoFeatures);
       
//       for (int site = 0; site < newBases.length; site+=defaultSiteCache) {
//            cacheTaxonSiteBlock(taxon, site);
//        }
       //myIsDirty=true;
    }
    
    @Override
    public void setDepthForAlleles(int taxon, int site, byte[] values) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setReferenceAllele(int site, byte diploidAllele) {
        myReference[site] = diploidAllele;
    }

    @Override
    public synchronized void addTaxon(Identifier id) {
        String basesPath = HapMapHDF5Constants.GENOTYPES + "/" + id.getFullName();
        myWriter.createByteArray(basesPath, myNumSites, genoFeatures);
        byte[] unkArray=new byte[getSiteCount()];
        Arrays.fill(unkArray, UNKNOWN_DIPLOID_ALLELE);
        myWriter.writeByteArray(basesPath, unkArray);
        myIdentifiers.add(id);
        myIdGroup=null;
        myIsDirty=true;
    }
    
    public synchronized void addTaxon(Identifier id, byte[] genotype, byte[][] depth) {
        int chunk=1<<16;
        if(myNumSites<chunk) chunk=myNumSites;
        String basesPath = HapMapHDF5Constants.GENOTYPES + "/" + id.getFullName();
        if(myWriter.exists(basesPath)) throw new IllegalStateException("Taxa Name Already Exists:"+basesPath);
        if(genotype.length!=myNumSites) throw new IllegalStateException("Setting all genotypes in addTaxon.  Wrong number of sites");
        
        myWriter.createByteArray(basesPath, myNumSites, chunk, genoFeatures);
//        myWriter.createByteArray(basesPath, myNumSites, genoFeatures);
//        myWriter.writeByteArray(basesPath, genotype, genoFeatures);
        int blocks=((myNumSites-1)/chunk)+1;
        for (int i = 0; i < blocks; i++) {
            int startPos=i*chunk;
            int length=(i<(blocks-1))?chunk:myNumSites-startPos;
            byte[] sub=new byte[length];
            System.arraycopy(genotype, startPos, sub, 0, length);
            myWriter.writeByteArrayBlockWithOffset(basesPath, sub, sub.length, (long)startPos);
        }
        int taxonIndex=myIdentifiers.size();
        myIdentifiers.add(id);
        myIdGroup=null;
        if(depth!=null) {
            if(depth.length!=6) throw new IllegalStateException("Just set A, C, G, T, -, + all at once");
            if(depth[0].length!=myNumSites) throw new IllegalStateException("Setting all depth in addTaxon.  Wrong number of sites");
            myWriter.writeByteMatrix(getTaxaDepthPath(taxonIndex), depth, genoFeatures);
        }
        myIsDirty=true;
    }

    @Override
    public synchronized void setTaxonName(int taxon, Identifier id) {
        if (taxon >= myIdentifiers.size()) {
            throw new IllegalStateException("MutableBitNucleotideAlignmentHDF5: setTaxonName: this taxa index does not exist: " + taxon);
        }
        String currentPath = getTaxaGenoPath(taxon);
        String newPath = HapMapHDF5Constants.GENOTYPES + "/" + id.getFullName();
        myWriter.move(currentPath, newPath);
        myIdentifiers.set(taxon, id);
        myIdGroup=null;
        myIsDirty=true;
    }

    @Override
    public synchronized void removeTaxon(int taxon) {
        String currentPath = getTaxaGenoPath(taxon);
     //   System.out.println(currentPath);
        myWriter.delete(currentPath);
        myIdGroup=null;
        myIdentifiers.remove(taxon);
       // System.out.println(getTaxaGenoPath(taxon));
        myGenoCache.invalidateAll();
        myIsDirty=true;
    }
    
    private void loadSiteDescriptorsToMemory() {
        myAlleleCnt=myWriter.readIntMatrix(HapMapHDF5Constants.ALLELE_CNT);
        myAlleleFreqOrder=myWriter.readByteMatrix(HapMapHDF5Constants.ALLELE_FREQ_ORD);
        maf=myWriter.readFloatArray(HapMapHDF5Constants.MAF);
        siteCov=myWriter.readFloatArray(HapMapHDF5Constants.SITECOV);
    }

    @Override
    public void clean() {
        long time=System.currentTimeMillis();
        System.out.print("Starting clean ...");
        recalcLocusCoordinates();
        HDF5AlignmentAnnotator faCalc=new HDF5AlignmentAnnotator(this, fileName, 
                HDF5AlignmentAnnotator.AnnoType.ALLELEFreq);
        faCalc.run();
        System.out.printf("finished in %dms%n",System.currentTimeMillis()-time);
        myIsDirty = false;
        myNumSites -= myNumSitesStagedToRemove;
        myNumSitesStagedToRemove = 0;
        //save the out cache to file
        //re-sort the taxa names
        //
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        return 1-maf[site];
    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        return maf[site];
    }

    @Override
    public byte getMajorAllele(int site) {
        return myAlleleFreqOrder[0][site]; 
    }

    @Override
    public byte getMinorAllele(int site) {
        return myAlleleFreqOrder[1][site];
    }
    
    public float getSiteCoverage(int site) {
        return siteCov[site];
    }

    @Override
    public boolean isDirty() {
        return myIsDirty;
    }

    private void setDirty() {
        myIsDirty = true;
    }

    private void setClean() {
        myIsDirty = false;
    }


    public int getDefaultCacheSize() {
        return defaultCacheSize;
    }

    public void setDefaultCacheSize(int defaultCacheSize) {
        this.defaultCacheSize = defaultCacheSize;
        initGenotypeCache();
    }

    public int getDefaultSiteCache() {
        return defaultSiteCache;
    }

    public void setDefaultSiteCache(int defaultSiteCache) {
        this.defaultSiteCache = defaultSiteCache;
        initGenotypeCache();
    }

    public void setCommonAlleles(int site, byte[] values) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    protected void setCalcAlleleFreq(int[][] af, byte[][] afOrder,
            float[] maf, float[] paf, float[] coverage, float[] hets) {
        myWriter.writeIntMatrix(HapMapHDF5Constants.ALLELE_CNT, af);
        myWriter.writeByteMatrix(HapMapHDF5Constants.ALLELE_FREQ_ORD, afOrder);
        myWriter.writeFloatArray(HapMapHDF5Constants.MAF, maf);
        myWriter.writeFloatArray(HapMapHDF5Constants.SITECOV, paf);
        myWriter.writeFloatArray(HapMapHDF5Constants.TAXACOV, coverage);
        myWriter.writeFloatArray(HapMapHDF5Constants.TAXAHET, hets);
    }

    @Override
    public void addSite(int site) {
        throw new UnsupportedOperationException("Will not be supported");
    }

    @Override
    public void removeSite(int site) {
        throw new UnsupportedOperationException("Will not be supported.");
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

    @Override
    public void setSNPID(int site, String name) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    
    
}
class LRUCache<K, V> extends LinkedHashMap<K, V> {
    private final int limit;
    public LRUCache(int limit) {
      super(16, 0.75f, true);
      this.limit = limit;
    }
    
    @Override
    protected boolean removeEldestEntry(Map.Entry<K,V> eldest) {
      return size() > limit;
    }
}
