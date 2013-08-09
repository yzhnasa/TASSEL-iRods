/*
 *  AnnotatedTaxon
 */

package net.maizegenetics.pal.taxa;

import com.google.common.collect.HashMultimap;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.site.GeneralAnnotation;

/**
 *
 * @author Ed Buckler
 */
public final class AnnotatedTaxon extends Identifier implements GeneralAnnotation {
      private final String myParent1;  //generally female
      private final String myParent2;  //generally male
      private final float myInbreedF; //inbreeding coefficient
      private final byte mySex;  //0=both, 1=female, 2=male
      private final String myPedigree;
    //Custom annotation are stored in the map
    private final HashMultimap<String, Object> myAnnoMap;

    @Override
    public Object[] getAnnotation(String annoName) {
        if(myAnnoMap==null) return null;
//        switch (annoName) {  //TODO: uncomment once in Java 7
//            case "myParent1":return myLocus;
//            case "myParent2":return myPosition;
//            case "myInbreedF":return myCM;
//            case "mySex":return myStrand;
//            case "pedigree":return mySNPID;
//        }
        return myAnnoMap.get(annoName).toArray();
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        try{return myAnnoMap.get(annoName).toArray(new String[0]);}
        catch(Exception e) {
            return null;
        }
    }

    @Override
    public Double[] getQuantAnnotation(String annoName) {
        try{
            Object[] o=myAnnoMap.get(annoName).toArray();
            if((o == null)||(!(o[0] instanceof Number))) return null;
            Double[] d=new Double[o.length];
            int i=0;
            for (Object o1 : o) {d[i++]=((Number)o1).doubleValue();}
            return d;
        }catch(Exception e) {
            return null;
        }
    }

    /**
     * A builder for creating immutable CorePosition instances.
     *<p> Example:
     * <pre>   {@code
     * AnnotatedTaxon cp= new AnnotatedTaxon.Builder("Z001E0001:Line:mays:Zea")
     *   .inbreedF(0.99)
     *   .parents("B73","B97")
     *   .pedigree("(B73xB97)S6I1")
     *   .build();}</pre>
     * <p>This would create an annotatedTaxon.
     */
    public static class Builder {
        // Required parameters
        private final Identifier myTaxon;
        // Optional parameters - initialized to default values
        private String parent1=null;  //generally female
        private String parent2=null;  //generally male
        private float inbreedF=Float.NaN;
        private byte sex=0;  //0=both, 1=female, 2=male
        private String pedigree=null;
        private HashMultimap<String, Object> myAnnoMap=null;

        /**
         * Constructor for Builder, requires a Taxon object
         * @param aTaxon taxon object
         */
        public Builder(Identifier aTaxon) {
            this.myTaxon = aTaxon;
        }
        /**
         * Constructor for Builder, requires a Taxon name
         * @param aTaxonName name of the taxon
         */
        public Builder(String aTaxonName) {
            this.myTaxon = new Identifier(aTaxonName);
        }
        /**Set sex: 0=both, 1=female, 2=male (default=0 Both)*/
        public Builder sex(byte val) {sex = val; return this;}
        /**Set inbreeding coefficient (default=Float.NaN)*/
        public Builder inbreedF(float val) {inbreedF = val; return this;}
        /**Set text definition of parents (default=null)*/
        public Builder parents(String mom, String dad) {parent1 = mom; parent2=dad; return this;}
        /**Set text definition of pedigree (default=null)*/
        public Builder pedigree(String val) {pedigree = val; return this;}

        /**Add non-standard annotation*/
        public Builder addAnno(String key, String value) {
            if(myAnnoMap==null) {
                myAnnoMap=HashMultimap.create(2,1);
            }
            myAnnoMap.put(key,value);
            return this;
        }
        /**Add non-standard annotation*/
        public Builder addAnno(String key, Number value) {
            if(myAnnoMap==null) {
                myAnnoMap=HashMultimap.create(2,1);
            }
            myAnnoMap.put(key,value);
            return this;
        }
        public AnnotatedTaxon build() {
            return new AnnotatedTaxon(this);
        }
    }

    private AnnotatedTaxon(Builder builder) {
        super(builder.myTaxon.getFullName());
        this.myParent1=builder.parent1;
        this.myParent2=builder.parent2;
        this.mySex=builder.sex;
        this.myPedigree=builder.pedigree;
        this.myInbreedF=builder.inbreedF;
        this.myAnnoMap=builder.myAnnoMap;
    }
}
