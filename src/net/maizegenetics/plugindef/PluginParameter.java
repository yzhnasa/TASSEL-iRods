package net.maizegenetics.plugindef;

import com.google.common.collect.Range;

/**
 * Defines the attributes of parameters to be used in the plugins
 *
 * @author Ed Buckler
 */
public class PluginParameter {  //TODO concept make this generic
    private final String guiName;
    private final String guiUnits;
    private final String cmdLineName;
    private final String cmdLineLongName;
    private final String description;
    private final Class type;
    private final Range range;
    private final Object value; //TODO make this of the type generic

    private PluginParameter(String guiName, String guiUnits, String cmdLineName,
                           String cmdLineLongName, String description, Class type, Range range, Object value) {
        this.guiName=guiName;
        this.guiUnits=guiUnits;
        this.cmdLineName=cmdLineName;
        this.cmdLineLongName=cmdLineLongName;
        this.description=description;
        this.type=type;
        this.range=range;
        this.value=value;
    }

    /**
     * Use these to change the value of an existing parameter, e.g. after a user changes the value.  Otherwise use
     * the Builder to create the parameter
     * @param oldParameter
     * @param newValue
     */
    public PluginParameter(PluginParameter oldParameter, Object newValue) {
        this.guiName=oldParameter.guiName;
        this.guiUnits=oldParameter.guiUnits;
        this.cmdLineName=oldParameter.cmdLineName;
        this.cmdLineLongName=oldParameter.cmdLineLongName;
        this.description=oldParameter.description;
        this.type=oldParameter.type;
        this.range=oldParameter.range;
        this.value=newValue;
    }

    public String guiName() {
        return guiName;
    }

    public String guiUnits() {
        return guiUnits;
    }

    public String cmdLineName() {
        return cmdLineName;
    }

    public String cmdLineLongName() {
        return cmdLineLongName;
    }

    public String description() {
        return description;
    }

    public Class type() {
        return type;
    }

    public Range range() {
        return range;
    }

    public Object value() {
        return value;
    }   //TODO change to generic return


    public static class Builder { //TODO make generic
        private String guiName;
        private String guiUnits="";
        private String cmdLineName;
        private String cmdLineLongName="";
        private String description="";
        private Class type;     //TODO remove as generic should know its type
        private Range range=null;
        private Object value;    //TODO make generic

        public Builder(String guiName, String cmdLineName, Class type, Object value) {
            //TODO perhaps simplify to  cmdLineName and value, guiName will be cmdLineName by default
            this.guiName=guiName;
            this.cmdLineName=cmdLineName;
            this.type=type;
            this.value=value;
        }

        public Builder guiUnits(String guiUnits) {
            this.guiUnits=guiUnits;
            return this;
        }

        public Builder cmdLineLongName(String cmdLineLongName) {
            this.cmdLineLongName=cmdLineLongName;
            return this;
        }

        public Builder description(String description) {
            this.description=description;
            return this;
        }

        public Builder range(Range range) {
            this.range=range;
            return this;
        }

        public PluginParameter build() {
            if(cmdLineLongName.isEmpty())  cmdLineLongName=cmdLineName;
            if(description.isEmpty())  description=guiName;
            return new PluginParameter(guiName, guiUnits, cmdLineName,
                     cmdLineLongName,  description, type, range, value);
        }
    }
}
