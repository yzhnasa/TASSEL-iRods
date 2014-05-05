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
    private boolean myMustBeChanged;
    private final Class<T> myClass;

    private PluginParameterTerry(String guiName, String guiUnits, String cmdLineName,
            String description, Range<T> range, T value, boolean required, Class<T> type) {
        myGuiName = guiName;
        myUnits = guiUnits;
        myCmdLineName = cmdLineName;
        myDescription = description;
        myRange = range;
        myValue = value;
        if ((myRange != null) && (!myRange.contains(myValue))) {
            throw new IllegalArgumentException("PluginParameterTerry: init: " + myCmdLineName + " value: " + value.toString() + " outside range: " + myRange.toString());
        }
        myRequired = required;
        myMustBeChanged = required;
        myClass = type;
        ;
    }

    /**
     * Use these to change the value of an existing parameter, e.g. after a user
     * changes the value. Otherwise use the Builder to create the parameter
     *
     * @param oldParameter
     * @param newValue
     */
    public PluginParameterTerry(PluginParameterTerry<T> oldParameter, T newValue) {
        this(oldParameter.myGuiName, oldParameter.myUnits, oldParameter.myCmdLineName,
                oldParameter.myDescription, oldParameter.myRange, newValue,
                oldParameter.myRequired, oldParameter.myClass);
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

    public boolean required() {
        return myRequired;
    }

    public Class<T> valueType() {
        return myClass;
    }

    public static class Builder<T extends Comparable<T>> {

        private String myGuiName;
        private String myUnits = "";
        private final String myCmdLineName;
        private String myDescription = "";
        private Range<T> myRange = null;
        private final T myValue;
        private boolean myIsRequired = false;
        private final Class<T> myClass;

        public Builder(Enum cmdLineName, T defaultValue, Class<T> type) {
            myCmdLineName = cmdLineName.toString();
            myValue = defaultValue;
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

        public Builder<T> required(boolean required) {
            myIsRequired = required;
            return this;
        }

        public Builder<T> guiName(String guiName) {
            myGuiName = guiName;
            return this;
        }

        public PluginParameterTerry<T> build() {
            if ((myGuiName == null) || (myGuiName.isEmpty())) {
                StringBuilder builder = new StringBuilder();
                builder.append(Character.toUpperCase(myCmdLineName.charAt(0)));
                for (int i = 1; i < myCmdLineName.length(); i++) {
                    char current = myCmdLineName.charAt(i);
                    if (Character.isUpperCase(current)) {
                        builder.append(" ");
                    }
                    builder.append(current);
                }
                myGuiName = builder.toString();
            }
            if (myDescription.isEmpty()) {
                myDescription = myGuiName;
            }
            return new PluginParameterTerry<>(myGuiName, myUnits, myCmdLineName,
                    myDescription, myRange, myValue, myIsRequired, myClass);
        }
    }
}
