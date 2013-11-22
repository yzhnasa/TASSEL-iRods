/*
 * IdGroupUtils
 */
package net.maizegenetics.taxa;

import java.util.ArrayList;
import java.util.TreeSet;

/**
 *
 * @author terry
 * @deprecated
 */
@Deprecated
public final class IdGroupUtils {

    private IdGroupUtils() {
        // To prevent instantiation
    }

    /**
     * @return true if <i>sub</i> IdGroup completely contained within <i>full</i>, false otherwise
     */
    public static boolean isContainedWithin(TaxaList sub, TaxaList full) {
        for (int i = 0; i < sub.getTaxaCount(); i++) {
            boolean found = false;
            Taxon subID = sub.get(i);
            for (int j = 0; j < full.getTaxaCount(); j++) {
                Taxon fullID = full.get(j);
                if (fullID.equals(subID)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                return false;
            }
        }
        return true;
    }

    /**
     * @return true if <i>id1</i> and <i>id2</i> share exactly the same identifiers (.equals() based, not reference base). The order is not important.
     */
    public static boolean isEqualIgnoringOrder(TaxaList id1, TaxaList id2) {
        return (isContainedWithin(id1, id2) && isContainedWithin(id2, id1));
    }

    /**
     * A convenience implementation of whichIdNumber that can be used by
     * IdGroup implementations
     * @return -1 if <i>s</i> not in <i>group</i>
     */
    public static int whichIdNumber(TaxaList group, String s) {
        for (int i = 0; i < group.getTaxaCount(); i++) {
            if (s.equals(group.get(i).getName())) {
                return i;
            }
        }
        return -1;
    }

    /**
     * @param group1 an IdGroup
     * @param group2 another IdGroup
     * @return the Ids in the intersection of groups 1 and 2, sorted in ascending order
     */
    public static TaxaList getCommonIds(TaxaList group1, TaxaList group2) {
        return getCommonIds(new TaxaList[]{group1, group2});
    }

    /**
     * Intersect joins the specified groups.
     *
     * @param groups groups to join.
     * @return The ids from the intersect join, sorted in ascending order
     */
    public static TaxaList getCommonIds(TaxaList[] groups) {

        if ((groups == null) || (groups.length == 0)) {
            return null;
        }

        TreeSet<Taxon> intersectIds = new TreeSet<Taxon>(groups[0]);
//        for (int x = 0; x < groups[0].getTaxaCount(); x++) {
//            intersectIds.add(groups[0].get(x));
//        }
        for (int i = 1; i < groups.length; i++) {
//            List temp = new ArrayList();
//            for (int j = 0; j < groups[i].getTaxaCount(); j++) {
//                temp.add(groups[i].get(j));
//            }
            intersectIds.retainAll(groups[i]);
        }

//        Taxon[] ids = new Taxon[intersectIds.size()];
//        intersectIds.toArray(ids);

        return new TaxaListBuilder().addAll(intersectIds).build();

    }

    /**
     * @param group1	an IdGroup
     * @param group2	another IdGroup
     * @return	the Ids in the union of groups 1 and 2, sorted in ascending order
     */
    public static TaxaList getAllIds(TaxaList group1, TaxaList group2) {
        return getAllIds(new TaxaList[]{group1, group2});
    }

    /**
     * Union joins the specified groups.
     *
     * @param groups groups to join.
     * @return The ids from the union join, sorted in ascending order
     * TODO move to Taxalist builder
     */
    public static TaxaList getAllIds(TaxaList[] groups) {

        if ((groups == null) || (groups.length == 0)) {
            return null;
        }

        TreeSet<Taxon> allIds = new TreeSet<Taxon>();
        for (int i = 0; i < groups.length; i++) {
            int n = groups[i].getTaxaCount();
            for (int j = 0; j < n; j++) {
                allIds.add(groups[i].get(j));
            }
        }
        return new TaxaListBuilder().addAll(allIds).build();

    }

    /**
     * @param original an IdGroup
     * @param include a boolean array with the same number of elements as the IdGroup original.
     * If the value of include is true for the ith element then that Taxon will be included
     * in the new subset. Otherwise, it will not.
     * @return a new IdGroup that is a subset of the original
     * @throws IllegalArgumentException if the number of Identifiers in original is not equal to the length of include.
     * TODO move to Taxalist builder
     */
    public static TaxaList idGroupSubset(TaxaList original, boolean[] include) {
        int nOld = original.getTaxaCount();
        if (nOld != include.length) {
            throw new IllegalArgumentException("Size of IdGroup and include array are different.");
        }
        ArrayList<Taxon> newIds = new ArrayList<Taxon>();
        for (int i = 0; i < nOld; i++) {
            if (include[i]) {
                newIds.add(original.get(i));
            }
        }
        return new TaxaListBuilder().addAll(newIds).build();
        //new SimpleIdGroup(newIds.toArray(new Taxon[newIds.size()]));
    }

    /**
     * Translates an array of identifiers into an array of strings
     */
    public static String[] getNames(Taxon[] ids) {
        String[] names = new String[ids.length];
        for (int i = 0; i < names.length; i++) {
            names[i] = ids[i].getName();
        }
        return names;
    }

    /**
     * Translates an array of identifiers into an array of strings, with optional removal of particular identifier
     * @param toIgnore the index of an idetifier to ignore, if <0 no element is ignored
     */
    public static String[] getNames(Taxon[] ids, int toIgnore) {
        if (toIgnore < 0 || toIgnore >= ids.length) {
            return getNames(ids);
        }
        String[] names = new String[ids.length - 1];
        int index = 0;
        for (int i = 0; i < names.length; i++) {
            if (i != toIgnore) {
                names[index] = ids[i].getName();
                index++;
            }
        }
        return names;
    }

    /**
     * Translates an an array of strings into an array of identifiers
     */
    public static Taxon[] getIdentifiers(String[] names) {
        Taxon[] ids = new Taxon[names.length];
        for (int i = 0; i < names.length; i++) {
            ids[i] = new Taxon(names[i]);
        }
        return ids;
    }

    /**
     * Translates an IdGroup into an array of identifiers
     */
    public static Taxon[] getIdentifiers(TaxaList idGroup) {
        Taxon[] ids = new Taxon[idGroup.getTaxaCount()];
        for (int i = 0; i < ids.length; i++) {
            ids[i] = idGroup.get(i);
        }
        return ids;
    }

    /**
     * Translates an IdGroup into an array of strings
     */
    public static String[] getNames(TaxaList ids) {
        String[] names = new String[ids.getTaxaCount()];
        for (int i = 0; i < names.length; i++) {
            names[i] = ids.get(i).getName();
        }
        return names;
    }

    /**
     * Translates an IDgroup into an array of strings, with optional removal of particular identifier
     * @param toIgnore the index of an idetifier to ignore, if <0 no element is ignored
     */
    public static String[] getNames(TaxaList ids, int toIgnore) {
        if (toIgnore < 0 || toIgnore >= ids.getTaxaCount()) {
            return getNames(ids);
        }
        int numberOfIDS = ids.getTaxaCount();
        String[] names = new String[numberOfIDS - 1];
        int index = 0;
        for (int i = 0; i < numberOfIDS; i++) {
            if (i != toIgnore) {
                names[index] = ids.get(i).getName();
                index++;
            }
        }
        return names;
    }
}
