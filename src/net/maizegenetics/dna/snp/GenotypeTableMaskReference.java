/*
 * GenotypeTableMaskReference
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.taxa.Taxon;

import java.awt.*;

/**
 *
 * @author terry
 */
public class GenotypeTableMaskReference extends AbstractGenotypeTableMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private final int myTaxonReference;
    private final GenotypeTable myAlignment;

    private GenotypeTableMaskReference(GenotypeTable align, int taxonReference, String name, Color color) {
        super(align, name, color, GenotypeTableMask.MaskType.reference);
        myTaxonReference = taxonReference;
        myAlignment = align;
    }

    public static GenotypeTableMaskReference getInstanceCompareReference(GenotypeTable align) {
        return getInstanceCompareReference(align, -1);
    }

    public static GenotypeTableMaskReference getInstanceCompareReference(GenotypeTable align, Taxon id) {
        int index = align.taxa().indexMatchingTaxon(id);
        if (index< 0) {
            throw new IllegalArgumentException("AlignmentMaskReference: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
    }

    public static GenotypeTableMaskReference getInstanceCompareReference(GenotypeTable align, String id) {
        int index = align.taxa().indexMatchingTaxon(id);
        if (index <0) {
            throw new IllegalArgumentException("AlignmentMaskReference: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
    }

    public static GenotypeTableMaskReference getInstanceCompareReference(GenotypeTable align, int index) {
        if ((index < -1) || (index >= align.numberOfTaxa())) {
            throw new IllegalArgumentException("AlignmentMaskReference: getInstanceCompareReference: unknown index: " + index);
        }
        String name;
        if (index == -1) {
            name = "Alignment Reference";
        } else {
            name = align.taxaName(index) + " Reference";
        }
        return new GenotypeTableMaskReference(align, index, name, getNextColor());
    }

    @Override
    public byte getMask(int taxon, int site) {
        if ((myTaxonReference == -1) && (GenotypeTableUtils.isEqualOrUnknown(myAlignment.genotype(taxon, site), myAlignment.referenceGenotype(site)))) {
            return 0;
        } else if (GenotypeTableUtils.isEqualOrUnknown(myAlignment.genotypeArray(taxon, site), myAlignment.genotypeArray(myTaxonReference, site))) {
            return 0;
        } else {
            return 1;
        }
    }
}
