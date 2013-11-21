/*
 * MergeAlignmentsPlugin
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.List;

/**
 *
 * @author terry
 */
public class MergeAlignmentsPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(MergeAlignmentsPlugin.class);

    public MergeAlignmentsPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        List<Datum> inputs = input.getDataOfType(Alignment.class);

        if ((inputs == null) || (inputs.size() < 2)) {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "Must select at least two alignments.");
            } else {
                myLogger.warn("performFunction: Must select at least two alignments.");
            }
            return null;
        }

        try {
            Alignment[] alignments = new Alignment[inputs.size()];
            for (int i = 0; i < inputs.size(); i++) {
                alignments[i] = (Alignment) ((Datum) inputs.get(i)).getData();
            }

            Alignment alignment = getInstance(alignments);
            DataSet result = new DataSet(new Datum("Merged Alignment", alignment, null), this);

            fireDataSetReturned(new PluginEvent(result, MergeAlignmentsPlugin.class));

            return result;
        } finally {
            fireProgress(100);
        }

    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = SeparatePlugin.class.getResource("images/Merge.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Merge Alignments";
    }

    @Override
    public String getToolTipText() {
        return "Merge Alignments";
    }

    //todo TAS-54 covers some of these topics, as we want general merging rule sets
    public static Alignment getInstance(Alignment[] alignments) {
       throw new UnsupportedOperationException("Logic for merging can be done much better now in TASSEL 5");
//        if ((alignments == null) || (alignments.length == 0)) {
//            return null;
//        }
//
//        String[][] resultEncodings = alignments[0].getAlleleEncodings();
//        if (resultEncodings.length != 1) {
//            throw new IllegalArgumentException("MutableSingleEncodeAlignment: getInstance: Alignments must have single allele encoding.");
//        }
//
//        String[][][] encodings = new String[alignments.length][][];
//        for (int i = 0; i < alignments.length; i++) {
//            encodings[i] = alignments[i].getAlleleEncodings();
//        }
//        if (!AlignmentUtils.areEncodingsEqual(encodings)) {
//            throw new IllegalArgumentException("MutableSingleEncodeAlignment: getInstance: Alignments must have same allele encoding.");
//        }
//
//        TreeSet<Taxon> taxa = new TreeSet<>();
//        List<String> siteNames = new ArrayList<>();
//        List<Integer> physicalPositions = new ArrayList<Integer>();
//        List<Chromosome> locusToLociIndex = new ArrayList<>();
//        List<Integer> locusIndices = new ArrayList<Integer>();
//
//        for (int i = 0; i < alignments.length; i++) {
//
//            TaxaList currentIds = alignments[i].getTaxaList();
//            for (int j = 0, n = currentIds.getIdCount(); j < n; j++) {
//                Taxon current = currentIds.getIdentifier(j);
//                if (taxa.contains(current)) {
//                    Taxon match = taxa.floor(current);
//                    Taxon merged = Identifier.getMergedInstance(match, current);
//                    taxa.remove(match);
//                    taxa.add(merged);
//                } else {
//                    taxa.add(current);
//                }
//            }
//
//            for (int s = 0, m = alignments[i].getSiteCount(); s < m; s++) {
//                String currentSiteName = alignments[i].getSNPID(s);
//                int currentPhysicalPos = alignments[i].getPositionInLocus(s);
//                Chromosome currentLocus = alignments[i].getLocus(s);
//                int index = siteNames.indexOf(currentSiteName);
//                if (index == -1) {
//                    siteNames.add(currentSiteName);
//                    physicalPositions.add(currentPhysicalPos);
//                    //int locusIndex = locusToLociIndex.indexOf(currentLocus);
//                    int locusIndex = -1;
//                    for (int li = 0; li < locusToLociIndex.size(); li++) {
//                        if (currentLocus.getChromosomeName().equals(locusToLociIndex.get(li).getChromosomeName())) {
//                            locusIndex = li;
//                            locusToLociIndex.set(li, Locus.getMergedInstance(currentLocus, locusToLociIndex.get(li)));
//                            break;
//                        }
//                    }
//                    if (locusIndex == -1) {
//                        locusIndices.add(locusToLociIndex.size());
//                        locusToLociIndex.add(currentLocus);
//                    } else {
//                        locusIndices.add(locusIndex);
//                    }
//                } else {
//                    if (i == 0) {
//                        throw new IllegalStateException("MutableSingleEncodeAlignment: getInstance: Duplicate site name in alignment: " + currentSiteName);
//                    } else {
//                        if (currentPhysicalPos != physicalPositions.get(index)) {
//                            throw new IllegalStateException("MutableSingleEncodeAlignment: getInstance: Physical Positions do not match for site name: " + currentSiteName);
//                        }
//                        int locusIndex = -1;
//                        for (int li = 0; li < locusToLociIndex.size(); li++) {
//                            if (currentLocus.getChromosomeName().equals(locusToLociIndex.get(li).getChromosomeName())) {
//                                locusIndex = li;
//                                break;
//                            }
//                        }
//                        if (locusIndices.get(index) != locusIndex) {
//                            throw new IllegalStateException("MutableSingleEncodeAlignment: getInstance: Loci do not match for site name: " + currentSiteName + " expecting: " + locusToLociIndex.get(locusIndices.get(index)) + " but doesn't match: " + currentLocus);
//                        }
//                    }
//                }
//            }
//
//        }
//        int[] variableSites = new int[physicalPositions.size()];
//        for (int i = 0, n = physicalPositions.size(); i < n; i++) {
//            variableSites[i] = physicalPositions.get(i);
//        }
//
//        int[] locusIndicesArray = new int[locusIndices.size()];
//        for (int i = 0, n = locusIndices.size(); i < n; i++) {
//            locusIndicesArray[i] = (int) locusIndices.get(i);
//        }
//
//        String[] siteNamesArray = new String[siteNames.size()];
//        siteNames.toArray(siteNamesArray);
//
//        List taxaList = new ArrayList<Identifier>(taxa);
//        MutableSingleEncodeAlignment result = null;
//
//        encodings = new String[2][][];
//        encodings[0] = NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;
//        encodings[1] = resultEncodings;
//        if (AlignmentUtils.areEncodingsEqual(encodings)) {
//            result = MutableNucleotideAlignment.getInstance(taxaList, variableSites, locusToLociIndex, locusIndicesArray, siteNamesArray);
//        } else {
//            result = getInstance(resultEncodings, taxaList, variableSites, locusToLociIndex, locusIndicesArray, siteNamesArray);
//        }
//        result.sortSitesByPhysicalPositionEmptyData();
//        result.setClean();
//
//        for (int i = 0; i < alignments.length; i++) {
//            myLogger.info("Merging Alignment: " + (i + 1) + " of " + alignments.length);
//            Alignment currentAlignment = alignments[i];
//            IdGroup ids = currentAlignment.getTaxaList();
//            int numSeqs = ids.getIdCount();
//            int[] taxaIndices = new int[numSeqs];
//            for (int t = 0; t < numSeqs; t++) {
//                taxaIndices[t] = taxaList.indexOf(ids.getIdentifier(t));
//            }
//            for (int s = 0, n = currentAlignment.getSiteCount(); s < n; s++) {
//                String siteName = currentAlignment.getSNPID(s);
//                int physicalPosition = currentAlignment.getPositionInLocus(s);
//                Locus locus = currentAlignment.getLocus(s);
//                for (int li = 0; li < locusToLociIndex.size(); li++) {
//                    if (locus.getChromosomeName().equals(locusToLociIndex.get(li).getChromosomeName())) {
//                        locus = locusToLociIndex.get(li);
//                        break;
//                    }
//                }
//                int site = result.getSiteOfPhysicalPosition(physicalPosition, locus, siteName);
//                if (site < 0) {
//                    throw new IllegalStateException("MutableSingleEncodeAlignment: getInstance: physical position: " + physicalPosition + " in locus: " + locus.getName() + " not found.");
//                }
//
//                for (int t = 0; t < numSeqs; t++) {
//                    result.setBase(taxaIndices[t], site, currentAlignment.getBase(t, s));
//                }
//            }
//        }
//
//        return result;

    }
}
