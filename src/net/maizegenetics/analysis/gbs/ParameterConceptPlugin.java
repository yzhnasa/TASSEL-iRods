/*
 * ParameterConceptPlugin
 */
package net.maizegenetics.analysis.gbs;

import java.awt.Frame;
import java.util.LinkedHashMap;
import java.util.Map;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.log4j.Logger;

import net.maizegenetics.plugindef.PluginParameterTerry;

/**
 *
 * @author Terry Casstevens
 */
public class ParameterConceptPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ParameterConceptPlugin.class);

    private final Map<PARAMETERS, PluginParameterTerry<?>> myParameters = new LinkedHashMap<>();

    public static enum PARAMETERS {

        inputFile;

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
        addStringParameter("Input File", PARAMETERS.inputFile, null, true);
    }

    private void addStringParameter(String guiName, PARAMETERS key, String defaultValue, boolean isRequired) {
        PluginParameterTerry<String> parameter = new PluginParameterTerry.Builder<>(guiName, key.toString(), defaultValue, isRequired).build();
        myParameters.put(key, parameter);
    }

    @Override
    public void setParameters(String[] args) {

        if (args == null || args.length == 0) {
            printUsage();
            System.exit(1);
        }

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            if (args[i].startsWith("-")) {
                arg = arg.substring(1);
                PARAMETERS key = PARAMETERS.valueOf(arg);
                PluginParameterTerry<?> parameter = myParameters.get(key);
                if ((i != args.length - 1) && (args[i + 1]).startsWith("-")) {
                    PluginParameterTerry newParameter = new PluginParameterTerry(parameter, true);
                } else {
                    PluginParameterTerry newParameter = new PluginParameterTerry(parameter, args[i + 1]);
                    myParameters.put(key, newParameter);
                    i++;
                }
            } else {
                throw new IllegalArgumentException("AbstractPlugin:  ");
            }
        }

    }

    private void printUsage() {
        // Get This from myParameters
        myLogger.info(
                "\nUsage:\n"
                + "ParameterConceptPluginPlugin <options>\n"
                + " -i  Input File Name\n");
    }

    public ParameterConceptPlugin setInput(String filename) {
        return this;
    }

    public String getInput() {
        return null;
    }

    public Object getParameterValue(PARAMETERS key) {
        return myParameters.get(key);
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
