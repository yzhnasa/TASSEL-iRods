
package net.maizegenetics.pal.ids;

import java.util.HashMap;


/**
 * 
 * @author Ed Buckler
 */

public class AnnotatedIdentifier extends Identifier {
    private HashMap<String, String> textAnnotations=null;
    private HashMap<String, Double> quantAnnotations=null;

    public AnnotatedIdentifier(String name) {
        super(name);
        
    }
    
    public AnnotatedIdentifier(String name, HashMap<String, String> textAnnotations, 
            HashMap<String, Double> quantAnnotations) {
        super(name);
        this.textAnnotations=textAnnotations;
        this.quantAnnotations=quantAnnotations;
    }
    
    public void addAnnotation(String annoName, String value) {
        textAnnotations.put(annoName, value);
    }
    
    public void addAnnotation(String annoName, double value) {
        quantAnnotations.put(annoName, value);
    }
    
    public String getTextAnnotation(String annoName) {
        return textAnnotations.get(annoName);
    }
    
    public double getQuantAnnotation(String annoName) {
        return quantAnnotations.get(annoName);
    }
    
    
}
