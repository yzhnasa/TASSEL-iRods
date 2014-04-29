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

    private final Map<String, PluginParameterTerry<?>> myParameters = new LinkedHashMap<>();

    public static enum PARAMETERS {

        inputFile, useReference, minAlleleFreq;

    };

    public ParameterConceptPlugin(Frame parentFrame) {
        super(parentFrame, false);
        defineParameters();
    }

    @Override
    public DataSet performFunction(DataSet input) {

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
                    if ((i != args.length - 1) && (args[i + 1]).startsWith("-")) {
                        PluginParameterTerry<Boolean> newParameter = new PluginParameterTerry(parameter, Boolean.TRUE);
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

    private void printUsage() {
        StringBuilder builder = new StringBuilder();
        builder.append("\nUsage:\n");
        builder.append(Utils.getBasename(getClass().getName())).append(" <options>\n");
        for (PluginParameterTerry<?> current : myParameters.values()) {
            builder.append("-");
            builder.append(current.cmdLineName());
            builder.append(" ");
            builder.append(current.description());
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
