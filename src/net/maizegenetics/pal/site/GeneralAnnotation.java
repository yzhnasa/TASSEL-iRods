package net.maizegenetics.pal.site;

/**
 * Provide generalized annotations (descriptors) for a taxon or site.
 *
 */
public interface GeneralAnnotation {

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (null if not present)
     */
    public Object[] getAnnotation(String annoName);

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (null if not present or numeric)
     */
    public String[] getTextAnnotation(String annoName);

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (null if not present or text)
     */
    public Double[] getQuantAnnotation(String annoName);

    //should we provide methods, to average the quantitative annotations, the first annotation
    //
}
