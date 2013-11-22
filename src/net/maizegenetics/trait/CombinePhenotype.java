package net.maizegenetics.trait;

import net.maizegenetics.trait.AbstractPhenotype;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import net.maizegenetics.pal.taxa.Taxon;

import java.util.*;
import java.util.Map.Entry;

public class CombinePhenotype extends AbstractPhenotype {

    Phenotype[] phenotypes;
    int numberOfPhenotypes;
    boolean isUnion;
    int[][] phenotypeColumn; //first dimension is columns for this Phenotype, second dimension is phenotype index then offset
    int[][] phenotypeRow;	//first dimension is rows for this Phenotype, second dimension are rows in sub-phenotypes
    int totalColumns;
    int totalRows;

    private CombinePhenotype(Phenotype[] phenotypes, boolean isUnion, TaxaList taxaGroup, List<Trait> traitList, int[][] rowmap, int[][] columnmap) {
        super(taxaGroup, traitList);
        this.phenotypes = phenotypes;
        this.isUnion = isUnion;
        phenotypeRow = rowmap;
        phenotypeColumn = columnmap;
        totalColumns = phenotypeColumn.length;
        totalRows = phenotypeRow.length;
    }

    public static CombinePhenotype getInstance(Phenotype phenotype1, Phenotype phenotype2, boolean isUnion) {
        return getInstance(new Phenotype[]{phenotype1, phenotype2}, isUnion);
    }

    public static CombinePhenotype getInstance(Phenotype[] phenotypes, boolean isUnion) {
        Object[] rowinfo = mapRows(phenotypes, isUnion);
        Object[] colinfo = mapColumns(phenotypes, isUnion);

        int[][] rowmap = (int[][]) rowinfo[0];
        TaxaList taxa = (TaxaList) rowinfo[1];
        int[][] colmap = (int[][]) colinfo[0];
        ArrayList<Trait> traitList = (ArrayList<Trait>) colinfo[1];

        return new CombinePhenotype(phenotypes, isUnion, taxa, traitList, rowmap, colmap);
    }

    private static Object[] mapRows(Phenotype[] phenotypes, boolean isUnion) {
        TreeMap<Taxon, int[]> rowTreeMap = new TreeMap<Taxon, int[]>();
        int npheno = phenotypes.length;
        int pcount = 0;
        for (Phenotype pheno : phenotypes) {
            int n = pheno.getNumberOfTaxa();
            for (int t = 0; t < n; t++) {
                Taxon id = pheno.getTaxon(t);
                int[] rows = rowTreeMap.get(id);
                if (rows == null) {
                    rows = new int[npheno];
                    for (int i = 0; i < npheno; i++) {
                        rows[i] = -1;
                    }
                    rowTreeMap.put(id, rows);
                }
                rows[pcount] = t;
            }
            pcount++;
        }

        //if intersect, delete taxa that are not in all phenotypes
        Set<Entry<Taxon, int[]>> rowSet = rowTreeMap.entrySet();
        if (!isUnion) {
            Iterator<Entry<Taxon, int[]>> rit = rowSet.iterator();
            while (rit.hasNext()) {
                boolean notComplete = false;
                int[] rows = rit.next().getValue();
                for (int i : rows) {
                    if (i == -1) {
                        notComplete = true;
                        break;
                    }
                }
                if (notComplete) {
                    rit.remove();
                }
            }
        }

        int nrows = rowSet.size();
        int[][] rowmap = new int[nrows][];
        Taxon[] ids = new Taxon[nrows];
        int count = 0;
        for (Entry<Taxon, int[]> entry : rowSet) {
            rowmap[count] = entry.getValue();
            ids[count] = entry.getKey();
            count++;
        }
        TaxaList tl=new TaxaListBuilder().addAll(ids).build();
        return new Object[]{rowmap, tl};
    }

    private static Object[] mapColumns(Phenotype[] phenotypes, boolean isUnion) {
        int ncol = 0;
        for (Phenotype pheno : phenotypes) {
            ncol += pheno.getNumberOfTraits();
        }
        ArrayList<Trait> traitList = new ArrayList<Trait>(ncol);
        int[][] colmap = new int[ncol][2];
        int count1 = 0;
        int count2 = 0;
        for (Phenotype pheno : phenotypes) {
            int n = pheno.getNumberOfTraits();
            for (int i = 0; i < n; i++) {
                colmap[count2][0] = count1;
                colmap[count2][1] = i;
                count2++;
                traitList.add(Trait.getInstance(pheno.getTrait(i)));
            }
            count1++;
        }

        return new Object[]{colmap, traitList};
    }

    //implement Phenotype methods
    public double[][] getData() {
        int nrows = getNumberOfTraits();
        int ncols = getNumberOfTraits();
        double[][] result = new double[nrows][ncols];
        for (int r = 0; r < nrows; r++) {
            for (int c = 0; c < ncols; c++) {
                result[r][c] = getData(r, c);
            }
        }
        return null;
    }

    public double getData(Taxon taxon, Trait trait) {
        return getData(whichTaxon(taxon), whichTrait(trait));
    }

    public double getData(int taxon, int trait) {
        int pheno = phenotypeColumn[trait][0];
        if (phenotypeRow[taxon][pheno] == -1) {
            return Double.NaN;
        } else {
            return phenotypes[pheno].getData(phenotypeRow[taxon][pheno], phenotypeColumn[trait][1]);
        }
    }

    public void setData(Taxon taxon, Trait trait, double value) {
        setData(whichTaxon(taxon), whichTrait(trait), value);
    }

    public void setData(int taxon, int trait, double value) {
        int pheno = phenotypeColumn[trait][0];
        phenotypes[pheno].setData(phenotypeRow[taxon][pheno], phenotypeColumn[trait][1], value);
    }

    public SimplePhenotype simpleCopy() {
        return new SimplePhenotype(getTaxa(), getTraits(), getData());
    }
}
