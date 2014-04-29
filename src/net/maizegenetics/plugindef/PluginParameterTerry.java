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
    private final boolean myRequired;
    private final boolean myMustBeChanged;
    private final Class<T> myClass;

    private PluginParameterTerry(String guiName, String guiUnits, String cmdLineName,
            String description, Range<T> range, T value, boolean required, Class<T> theClass) {
        myGuiName = guiName;
        myUnits = guiUnits;
        myCmdLineName = cmdLineName;
        myDescription = description;
        myRange = range;
        myValue = value;
        if (!myRange.contains(myValue)) {
            throw new IllegalArgumentException("PluginParameterTerry: init: new value: " + value.toString() + " outside range: " + myRange.toString());
        }
        myRequired = required;
        myMustBeChanged = required;
        myClass = theClass;
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
        if (!myRange.contains(myValue)) {
            throw new IllegalArgumentException("PluginParameterTerry: init: new value: " + newValue.toString() + " outside range: " + myRange.toString());
        }
        myRequired = oldParameter.myRequired;
        myMustBeChanged = false;
        myClass = oldParameter.myClass;
    }

    public PluginParameterTerry(PluginParameterTerry<T> oldParameter, String newValue) {
        myGuiName = oldParameter.myGuiName;
        myUnits = oldParameter.myUnits;
        myCmdLineName = oldParameter.myCmdLineName;
        myDescription = oldParameter.myDescription;
        myRange = oldParameter.myRange;
        myClass = oldParameter.myClass;
        myValue = convert(newValue, myClass);
        if (!myRange.contains(myValue)) {
            throw new IllegalArgumentException("PluginParameterTerry: init: new value: " + newValue + " outside range: " + myRange.toString());
        }
        myRequired = oldParameter.myRequired;
        myMustBeChanged = false;
    }

    private T convert(String input, Class<T> outputClass) {

        try {
            return input == null ? null : outputClass.getConstructor(String.class).newInstance(input);
        } catch (Exception e) {
            throw new IllegalArgumentException("PluginParameterTerry: convert: Unknown type: " + outputClass.getName());
        }

        // We might need something like this if above
        // doesn't handle all cases.
        // if (outputClass.isAssignableFrom(Double.class))
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

    public boolean required() {
        return myRequired;
    }

    public static class Builder<T extends Comparable<T>> {

        private final String myGuiName;
        private String myUnits = "";
        private final String myCmdLineName;
        private String myDescription = "";
        private Range<T> myRange = null;
        private final T myValue;
        private final boolean myIsRequired;
        private final Class<T> myClass;

        public Builder(String guiName, String cmdLineName, T value, boolean isRequired, Class<T> type) {
            myGuiName = guiName;
            myCmdLineName = cmdLineName;
            myValue = value;
            myIsRequired = isRequired;
            myClass = type;
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
            return new PluginParameterTerry<>(myGuiName, myUnits, myCmdLineName,
                    myDescription, myRange, myValue, myIsRequired, myClass);
        }
    }
}
