package net.maizegenetics.pal.ids;

import java.util.HashMap;
import java.io.Serializable;


/**
 * Created by IntelliJ IDEA.
 * User: PeterLocal
 * Date: Jan 30, 2007
 * Time: 11:30:04 AM
 *
 */
public class TraitId implements Serializable {
    protected HashMap properties;

    public TraitId(String name) {
        properties = new HashMap();
        properties.put("Trait", name);
    }

    public TraitId(String name, String environmentName) {
        properties = new HashMap();
        properties.put("Trait", name);
        properties.put("Env", environmentName);
    }

    public TraitId() {
        properties = new HashMap();
    }

    public TraitId(HashMap properties) {
        this.properties = properties;
    }

    public boolean equals(Object obj) {
        if(obj instanceof TraitId) {
            TraitId aTrait = (TraitId) obj;
            if (properties.equals(aTrait.properties)) return true;
        }

        return false;
    }

    public int hashCode() {
        return properties.hashCode();
    }

    public String toString() {
        return getName();
    }

    public String getName() {
        return (String) properties.get("Trait");
    }

    public String getEnvironmentName() {
        Object obj = properties.get("Env");
        if (obj instanceof String) return (String) obj;
        return "NA";
    }

    public Object getProperty(Object key) {
        return properties.get(key);
    }

    public void setProperty(Object key, Object value) {
        properties.put(key, value);
    }

    public HashMap getPropertyMap() {
        return properties;
    }

    public TraitId copy() {
        return new TraitId(new HashMap(properties));
    }

    public void changeKey(Object from, Object to) {
        if (properties.containsKey(from)) {
            properties.put(to, properties.get(from));
            properties.remove(from);
        }
    }

    public boolean containsKey(Object key) {
        return properties.containsKey(key);
    }
}
