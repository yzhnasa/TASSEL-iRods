package net.maizegenetics.plugindef;

import com.google.common.collect.Range;

/**
 * Defines the attributes of parameters to be used in the plugins
 *
 * @param <T>
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 *
 */
public class PluginParameterTerry<T extends Comparable<T>> {

    private final String myGuiName;
    private final String myUnits;
    private final String myCmdLineName;
    private final String myDescription;
    private final Range<T> myRange;
    private final T myValue;
    private final boolean myMustBeChanged;

    private PluginParameterTerry(String guiName, String guiUnits, String cmdLineName,
            String description, Range<T> range, T value, boolean mustBeChanged) {
        myGuiName = guiName;
        myUnits = guiUnits;
        myCmdLineName = cmdLineName;
        myDescription = description;
        myRange = range;
        myValue = value;
        myMustBeChanged = mustBeChanged;
    }

    /**
     * Use these to change the value of an existing parameter, e.g. after a user
     * changes the value. Otherwise use the Builder to create the parameter
     *
     * @param oldParameter
     * @param newValue
     */
    public PluginParameterTerry(PluginParameterTerry<T> oldParameter, T newValue) {
        myGuiName = oldParameter.myGuiName;
        myUnits = oldParameter.myUnits;
        myCmdLineName = oldParameter.myCmdLineName;
        myDescription = oldParameter.myDescription;
        myRange = oldParameter.myRange;
        myValue = newValue;
        myMustBeChanged = false;
    }
    
    public String guiName() {
        return myGuiName;
    }

    public String units() {
        return myUnits;
    }

    public String cmdLineName() {
        return myCmdLineName;
    }

    public String description() {
        return myDescription;
    }

    public Range<T> range() {
        return myRange;
    }

    public T value() {
        return myValue;
    }

    public boolean mustBeChanged() {
        return myMustBeChanged;
    }

    public static class Builder<T extends Comparable<T>> {

        private final String myGuiName;
        private String myUnits = "";
        private final String myCmdLineName;
        private String myDescription = "";
        private Range<T> myRange = null;
        private final T myValue;
        private final boolean myIsRequired;

        public Builder(String guiName, String cmdLineName, T value, boolean isRequired) {
            myGuiName = guiName;
            myCmdLineName = cmdLineName;
            myValue = value;
            myIsRequired = isRequired;
        }

        public Builder<T> units(String units) {
            myUnits = units;
            return this;
        }

        public Builder<T> description(String description) {
            myDescription = description;
            return this;
        }

        public Builder<T> range(Range<T> range) {
            myRange = range;
            return this;
        }

        public PluginParameterTerry<T> build() {
            if (myDescription.isEmpty()) {
                myDescription = myGuiName;
            }
            return new PluginParameterTerry<T>(myGuiName, myUnits, myCmdLineName,
                    myDescription, myRange, myValue, myIsRequired);
        }
    }
}
