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
public class PluginParameter<T extends Comparable<T>> {

    private final String myGuiName;
    private final String myUnits;
    private final String myCmdLineName;
    private final String myDescription;
    private final Range<T> myRange;
    private final T myDefaultValue;
    private final T myValue;
    private final boolean myRequired;
    private final Class<T> myClass;

    public enum FILE_TYPE {

        NA, IN, OUT
    };
    private final FILE_TYPE myFileType;

    private PluginParameter(String guiName, String guiUnits, String cmdLineName,
            String description, Range<T> range, T defaultValue, T value, boolean required, FILE_TYPE fileType, Class<T> type) {
        myGuiName = guiName;
        myUnits = guiUnits;
        myCmdLineName = cmdLineName;
        myDescription = description;
        myRange = range;
        myDefaultValue = defaultValue;
        if (value == null) {
            myValue = defaultValue;
        } else {
            myValue = value;
        }
        if ((myRange != null) && (!myRange.contains(myValue))) {
            throw new IllegalArgumentException("PluginParameter: init: " + myCmdLineName + " value: " + value.toString() + " outside range: " + myRange.toString());
        }
        myRequired = required;
        if ((myDefaultValue != null) && (myRequired)) {
            throw new IllegalArgumentException("PluginParameter: init: " + myCmdLineName + " shouldn't have default value and be required.");
        }
        myClass = type;
        myFileType = fileType;
    }

    /**
     * Use these to change the value of an existing parameter, e.g. after a user
     * changes the value. Otherwise use the Builder to create the parameter
     *
     * @param oldParameter
     * @param newValue
     */
    public PluginParameter(PluginParameter<T> oldParameter, T newValue) {
        this(oldParameter.myGuiName, oldParameter.myUnits, oldParameter.myCmdLineName,
                oldParameter.myDescription, oldParameter.myRange, oldParameter.myDefaultValue, newValue,
                oldParameter.myRequired, oldParameter.myFileType, oldParameter.myClass);
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

    public T defaultValue() {
        return myDefaultValue;
    }

    public boolean required() {
        return myRequired;
    }

    public Class<T> valueType() {
        return myClass;
    }

    public FILE_TYPE fileType() {
        return myFileType;
    }

    public static class Builder<T extends Comparable<T>> {

        private String myGuiName;
        private String myUnits = "";
        private final String myCmdLineName;
        private String myDescription = "";
        private Range<T> myRange = null;
        private final T myDefaultValue;
        private boolean myIsRequired = false;
        private final Class<T> myClass;
        private FILE_TYPE myFileType = FILE_TYPE.NA;

        public Builder(Enum cmdLineName, T defaultValue, Class<T> type) {
            myCmdLineName = cmdLineName.toString();
            myDefaultValue = defaultValue;
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

        public Builder<T> inFile() {
            myFileType = FILE_TYPE.IN;
            return this;
        }

        public Builder<T> outFile() {
            myFileType = FILE_TYPE.OUT;
            return this;
        }

        public PluginParameter<T> build() {
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
            return new PluginParameter<>(myGuiName, myUnits, myCmdLineName,
                    myDescription, myRange, myDefaultValue, null, myIsRequired, myFileType, myClass);
        }
    }
}
