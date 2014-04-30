/*
 * ParameterConceptPlugin
 */
package net.maizegenetics.analysis.gbs;

import com.google.common.collect.Range;

import java.awt.Frame;
import javax.swing.ImageIcon;

import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameterTerry;
import net.maizegenetics.util.Utils;

/**
 *
 * @author Terry Casstevens
 */
public class ParameterConceptPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ParameterConceptPlugin.class);

    private final Map<String, PluginParameterTerry<? extends Comparable>> myParameters = new LinkedHashMap<>();

    public static enum PARAMETERS {

        inputFile, useReference, minAlleleFreq;

    };

    public ParameterConceptPlugin(Frame parentFrame) {
        super(parentFrame, false);
        defineParameters();
    }

    @Override
    public DataSet performFunction(DataSet input) {

        printParameterValues();

        try {
            return null;
        } finally {
            fireProgress(100);
        }

    }

    private void defineParameters() {
        addStringParameter("Input File", PARAMETERS.inputFile.toString(), null, true, null, null, null);
        addBooleanParameter("Use Reference", PARAMETERS.useReference.toString(), Boolean.FALSE, false, null, null);
        addDoubleParameter("Minimum Allele Frequency", PARAMETERS.minAlleleFreq.toString(), 0.01, false, null, Range.closed(0.0, 1.0), "Percentage");
    }

    protected void addStringParameter(String guiName, String cmdLineName, String defaultValue, boolean isRequired, String description, Range<String> range, String units) {
        PluginParameterTerry.Builder<String> builder = new PluginParameterTerry.Builder<>(guiName, cmdLineName, defaultValue, isRequired, String.class);
        modifyBuilder(builder, description, range, units);
        myParameters.put(cmdLineName, builder.build());
    }

    protected void addBooleanParameter(String guiName, String cmdLineName, Boolean defaultValue, boolean isRequired, String description, String units) {
        PluginParameterTerry.Builder<Boolean> builder = new PluginParameterTerry.Builder<>(guiName, cmdLineName, defaultValue, isRequired, Boolean.class);
        modifyBuilder(builder, description, null, units);
        myParameters.put(cmdLineName, builder.build());
    }

    protected void addDoubleParameter(String guiName, String cmdLineName, Double defaultValue, boolean isRequired, String description, Range<Double> range, String units) {
        PluginParameterTerry.Builder<Double> builder = new PluginParameterTerry.Builder<>(guiName, cmdLineName, defaultValue, isRequired, Double.class);
        modifyBuilder(builder, description, range, units);
        myParameters.put(cmdLineName, builder.build());
    }
    
    protected void addFloatParameter(String guiName, String cmdLineName, Float defaultValue, boolean isRequired, String description, Range<Float> range, String units) {
        PluginParameterTerry.Builder<Float> builder = new PluginParameterTerry.Builder<>(guiName, cmdLineName, defaultValue, isRequired, Float.class);
        modifyBuilder(builder, description, range, units);
        myParameters.put(cmdLineName, builder.build());
    }
    
    protected void addIntegerParameter(String guiName, String cmdLineName, Integer defaultValue, boolean isRequired, String description, Range<Integer> range, String units) {
        PluginParameterTerry.Builder<Integer> builder = new PluginParameterTerry.Builder<>(guiName, cmdLineName, defaultValue, isRequired, Integer.class);
        modifyBuilder(builder, description, range, units);
        myParameters.put(cmdLineName, builder.build());
    }
    
    protected void addByteParameter(String guiName, String cmdLineName, Byte defaultValue, boolean isRequired, String description, Range<Byte> range, String units) {
        PluginParameterTerry.Builder<Byte> builder = new PluginParameterTerry.Builder<>(guiName, cmdLineName, defaultValue, isRequired, Byte.class);
        modifyBuilder(builder, description, range, units);
        myParameters.put(cmdLineName, builder.build());
    }
    
    protected void addLongParameter(String guiName, String cmdLineName, Long defaultValue, boolean isRequired, String description, Range<Long> range, String units) {
        PluginParameterTerry.Builder<Long> builder = new PluginParameterTerry.Builder<>(guiName, cmdLineName, defaultValue, isRequired, Long.class);
        modifyBuilder(builder, description, range, units);
        myParameters.put(cmdLineName, builder.build());
    }
    
    protected void addCharParameter(String guiName, String cmdLineName, Character defaultValue, boolean isRequired, String description, Range<Character> range, String units) {
        PluginParameterTerry.Builder<Character> builder = new PluginParameterTerry.Builder<>(guiName, cmdLineName, defaultValue, isRequired, Character.class);
        modifyBuilder(builder, description, range, units);
        myParameters.put(cmdLineName, builder.build());
    }

    private static <T extends Comparable<T>> void modifyBuilder(PluginParameterTerry.Builder<T> builder, String description, Range<T> range, String units) {
        if (description != null) {
            builder.description(description);
        }
        if (range != null) {
            builder.range(range);
        }
        if (units != null) {
            builder.units(units);
        }
    }

    @Override
    public void setParameters(String[] args) {

        if (args != null) {

            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (args[i].startsWith("-")) {
                    arg = arg.substring(1);
                    PluginParameterTerry<?> parameter = myParameters.get(arg);
                    if (parameter == null) {
                        myLogger.error("Unrecognized argument: " + args[i]);
                        printUsage();
                        System.exit(1);
                    }
                    if ((i == args.length - 1) || (args[i + 1]).startsWith("-")) {
                        if (parameter.valueType().isAssignableFrom(Boolean.class)) {
                            PluginParameterTerry<Boolean> newParameter = new PluginParameterTerry(parameter, Boolean.TRUE);
                            myParameters.put(arg, newParameter);
                        } else {
                            myLogger.error("Parameter requires a value: " + args[i]);
                            printUsage();
                            System.exit(1);
                        }
                    } else {
                        PluginParameterTerry<?> newParameter = new PluginParameterTerry(parameter, args[i + 1]);
                        myParameters.put(arg, newParameter);
                        i++;
                    }
                } else {
                    myLogger.error("Argument expected to start with dash(-): " + args[i]);
                    printUsage();
                    System.exit(1);
                }
            }

        }

        for (PluginParameterTerry<?> current : myParameters.values()) {
            if (current.mustBeChanged()) {
                myLogger.error("flag -" + current.cmdLineName() + " is required.\n");
                printUsage();
                System.exit(1);
            }
        }

    }

    private void printParameterValues() {
        StringBuilder builder = new StringBuilder();
        builder.append("\n");
        builder.append(Utils.getBasename(getClass().getName()));
        builder.append(" Parameters\n");
        for (PluginParameterTerry<?> current : myParameters.values()) {
            builder.append(current.cmdLineName());
            builder.append(": ");
            builder.append(current.value());
            builder.append("\n");
        }
        myLogger.info(builder.toString());
    }

    private void printUsage() {
        StringBuilder builder = new StringBuilder();
        builder.append("\nUsage:\n");
        builder.append(Utils.getBasename(getClass().getName())).append(" <options>\n");
        for (PluginParameterTerry<?> current : myParameters.values()) {
            builder.append("-");
            builder.append(current.cmdLineName());
            builder.append(" ");
            if (current.valueType().isAssignableFrom(Boolean.class)) {
                builder.append("<true | false>");
            } else {
                builder.append("<");
                builder.append(current.guiName());
                builder.append(">");
            }
            builder.append(" : ");
            builder.append(current.description());
            if (current.range() != null) {
                builder.append(" ");
                builder.append(current.range().toString());
            }
            if (current.required()) {
                builder.append(" (required)");
            }
            builder.append("\n");
        }
        myLogger.info(builder.toString());
    }

    public ParameterConceptPlugin setInput(String filename) {
        return this;
    }

    public String getInput() {
        return null;
    }

    public Object getParameterValue(PARAMETERS key) {
        return myParameters.get(key.toString());
    }

    public Object getParameterValue(String key) {
        return myParameters.get(key);
    }

    public <T extends Comparable<T>> ParameterConceptPlugin setParameterValue(String key, T value) {
        PluginParameterTerry parameter = myParameters.get(key);
        if (parameter == null) {
            throw new IllegalArgumentException("AbstractPlugin: setParameterValue: Unknown Parameter: " + key);
        }
        PluginParameterTerry<T> newParameter = new PluginParameterTerry<>(parameter, value);
        myParameters.put(key, newParameter);
        return this;
    }
    
    public <T extends Comparable<T>> ParameterConceptPlugin setParameterValue(Enum key, T value) {
        PluginParameterTerry parameter = myParameters.get(key.toString());
        if (parameter == null) {
            throw new IllegalArgumentException("AbstractPlugin: setParameterValue: Unknown Parameter: " + key.toString());
        }
        PluginParameterTerry<T> newParameter = new PluginParameterTerry<>(parameter, value);
        myParameters.put(key.toString(), newParameter);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}