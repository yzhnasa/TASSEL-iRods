/*
 * ParameterConceptPlugin
 */
package net.maizegenetics.analysis.gbs;

import com.google.common.collect.Range;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

import javax.swing.ImageIcon;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.JOptionPane;
import javax.swing.JFileChooser;
import javax.swing.JCheckBox;
import javax.swing.JComponent;

import java.io.File;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameterTerry;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.prefs.TasselPrefs;

/**
 *
 * @author Terry Casstevens
 */
public class ParameterConceptPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ParameterConceptPlugin.class);

    public enum PARAMETERS {

        inputFile, useReference, minAlleleFreq
    };

    private PluginParameterTerry<String> myInputFile = new PluginParameterTerry.Builder<String>(PARAMETERS.inputFile, null, String.class).required(true).inFile().build();
    private PluginParameterTerry<Double> myMinAlleleFreq = new PluginParameterTerry.Builder<Double>(PARAMETERS.minAlleleFreq, 0.01, Double.class).range(Range.closed(0.0, 1.0)).units("Ratio").build();
    private PluginParameterTerry<Boolean> myUseReference = new PluginParameterTerry.Builder<Boolean>(PARAMETERS.useReference, false, Boolean.class).build();

    public ParameterConceptPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public ParameterConceptPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public static void main(String[] args) {

        java.util.Properties props = new java.util.Properties();
        props.setProperty("log4j.logger.net.maizegenetics", "INFO, stdout");
        props.setProperty("log4j.appender.stdout",
                "org.apache.log4j.ConsoleAppender");
        props.setProperty("log4j.appender.stdout.layout",
                "org.apache.log4j.TTCCLayout");
        PropertyConfigurator.configure(props);

        ParameterConceptPlugin plugin = new ParameterConceptPlugin.Builder(null)
                .setParameter(PARAMETERS.inputFile, "terry.txt")
                .setParameter(PARAMETERS.minAlleleFreq, 0.2)
                .setParameter(PARAMETERS.useReference, true).build();

        plugin.performFunction(null);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        try {

            if (isInteractive()) {
                if (!setParametersViaGUI()) {
                    return null;
                }
            }

            printParameterValues();
            checkParameters();

            DataSet output = processData(input);
            if (output != null) {
                fireDataSetReturned(new PluginEvent(output, getClass()));
            }
            return output;

        } catch (Exception e) {

            if (isInteractive()) {
                StringBuilder builder = new StringBuilder();
                builder.append(e.getMessage());
                builder.append("\n");
                String str = builder.toString();
                DialogUtils.showError(str, getParentFrame());
            } else {
                myLogger.error(e.getMessage());
                e.printStackTrace();
            }
            return null;

        } finally {
            fireProgress(100);
        }

    }

    private DataSet processData(DataSet input) {
        return null;
    }

    private List<PluginParameterTerry<?>> getParameterInstances() {

        List<PluginParameterTerry<?>> result = new ArrayList<>();
        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameterTerry.class)) {
                try {
                    PluginParameterTerry<?> parameter = (PluginParameterTerry) current.get(this);
                    result.add(parameter);
                } catch (Exception e) {
                    e.printStackTrace();
                    throw new IllegalArgumentException("AbstractPlugin: getParameterInstances: problem getting parameter instances");
                }

            }
        }

        return result;
    }

    private Field getParameterField(String key) {

        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameterTerry.class)) {
                try {
                    PluginParameterTerry<?> parameter = (PluginParameterTerry) current.get(this);
                    if (parameter.cmdLineName().equals(key)) {
                        return current;
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    throw new IllegalArgumentException("AbstractPlugin: getParameterField: problem with key: " + key);
                }

            }
        }

        throw new IllegalArgumentException("AbstractPlugin: getParameterField: unknown key: " + key);
    }

    private PluginParameterTerry<?> getParameterInstance(String key) {
        try {
            Field field = getParameterField(key);
            return (PluginParameterTerry) field.get(this);
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("AbstractPlugin: getParameterInstance: problem with key: " + key);
        }
    }

    private <T extends Comparable<T>> Class<T> getParameterGeneric(String key) {

        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameterTerry.class)) {
                try {
                    PluginParameterTerry<T> parameter = (PluginParameterTerry<T>) current.get(this);
                    if (parameter.cmdLineName().equals(key)) {
                        System.out.println("generic type: " + current.getGenericType());
                        return (Class<T>) current.getGenericType().getClass();
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    throw new IllegalArgumentException("AbstractPlugin: getParameter: problem with key: " + key);
                }

            }
        }

        throw new IllegalArgumentException("AbstractPlugin: getParameter: unknown key: " + key);
    }

    public static <T> T convert(String input, Class<T> outputClass) {
        try {
            return input == null ? null : outputClass.getConstructor(String.class).newInstance(input);
        } catch (InvocationTargetException nfe) {
            throw new IllegalArgumentException("convert: Problem converting: " + input + " to " + outputClass.getName());
        } catch (Exception e) {
            throw new IllegalArgumentException("convert: Unknown type: " + outputClass.getName());
        }
    }

    @Override
    public void setParameters(String[] args) {

        if (args != null) {

            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (args[i].startsWith("-")) {
                    arg = arg.substring(1);
                    PluginParameterTerry<?> parameter = getParameterInstance(arg);
                    if (parameter == null) {
                        myLogger.error("Unrecognized argument: " + args[i]);
                        printUsage();
                        System.exit(1);
                    }
                    if ((i == args.length - 1) || (args[i + 1]).startsWith("-")) {
                        if (parameter.valueType().isAssignableFrom(Boolean.class)) {
                            setParameter(arg, Boolean.TRUE);
                        } else {
                            myLogger.error("Parameter requires a value: " + args[i]);
                            printUsage();
                            System.exit(1);
                        }
                    } else {
                        setParameter(arg, args[i + 1]);
                        i++;
                    }
                } else {
                    myLogger.error("Argument expected to start with dash(-): " + args[i]);
                    printUsage();
                    System.exit(1);
                }
            }

        }

    }

    /**
     * Verification checks of parameters.
     */
    private void checkParameters() {

        for (PluginParameterTerry<?> current : getParameterInstances()) {

            if (current.required()) {
                Object value = current.value();
                if ((value == null) || (value.toString().trim().length() == 0)) {
                    if (isInteractive()) {
                        throw new IllegalStateException(current.guiName() + " must be defined.");
                    } else {
                        myLogger.error("-" + current.cmdLineName() + " is required.\n");
                        printUsage();
                        System.exit(1);
                    }
                }
            }

            if (current.fileType() == PluginParameterTerry.FILE_TYPE.IN) {
                String filename = current.value().toString();
                if (!new File(filename).exists()) {
                    if (isInteractive()) {
                        throw new IllegalStateException(current.guiName() + ": " + filename + " doesn't exist.");
                    } else {
                        myLogger.error("-" + current.cmdLineName() + ": " + filename + " doesn't exist\n");
                        printUsage();
                        System.exit(1);
                    }
                }
            }

            if (current.fileType() == PluginParameterTerry.FILE_TYPE.OUT) {
                String filename = current.value().toString();
                String outFolder = Utils.getDirectory(filename);
                File outDir = new File(outFolder);
                if (!outDir.isDirectory()) {
                    if (isInteractive()) {
                        throw new IllegalStateException(current.guiName() + ": Output Directory: " + outFolder + " doesn't exist.");
                    } else {
                        myLogger.error("-" + current.cmdLineName() + ": Output Directory: " + outFolder + " doesn't exist\n");
                        printUsage();
                        System.exit(1);
                    }
                }
            }

        }

    }

    private void printParameterValues() {
        StringBuilder builder = new StringBuilder();
        builder.append("\n");
        builder.append(Utils.getBasename(getClass().getName()));
        builder.append(" Parameters\n");
        for (PluginParameterTerry<?> current : getParameterInstances()) {
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
        for (PluginParameterTerry<?> current : getParameterInstances()) {
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

    public String inputFile() {
        return myInputFile.value();
    }

    public Object getParameter(Enum key) {
        return getParameterInstance(key.toString()).value();
    }

    protected <T extends Comparable<T>> ParameterConceptPlugin setParameter(String key, T value) {

        PluginParameterTerry<T> parameter = null;
        try {

            Field field = getParameterField(key);
            parameter = (PluginParameterTerry<T>) field.get(this);
            if (parameter == null) {
                throw new IllegalArgumentException("setParameter: Unknown Parameter: " + key);
            }
            if ((parameter.range() != null) && (!parameter.range().contains(value))) {
                throw new IllegalArgumentException("setParameter: " + parameter.cmdLineName() + " value: " + value.toString() + " outside range: " + parameter.range().toString());
            }
            PluginParameterTerry<T> newParameter = new PluginParameterTerry<>(parameter, value);
            field.set(this, newParameter);

        } catch (Exception e) {
            if (isInteractive()) {
                try {
                    throw e;
                } catch (IllegalAccessException ex) {
                    ex.printStackTrace();
                }
            } else {
                myLogger.error(key + ": " + e.getMessage());
                printUsage();
                System.exit(1);
            }
        }

        return this;
    }

    protected <T extends Comparable<T>> ParameterConceptPlugin setParameter(String key, String value) {

        PluginParameterTerry<T> parameter = null;
        try {
            System.out.println("generic: " + getParameterGeneric(key));
            parameter = (PluginParameterTerry<T>) getParameterInstance(key);
            return setParameter(key, convert(value, parameter.valueType()));
        } catch (Exception e) {
            if (isInteractive()) {
                throw e;
            } else {
                myLogger.error(key + ": " + e.getMessage());
                printUsage();
                System.exit(1);
            }
        }
        return this;
    }

    protected <T extends Comparable<T>> ParameterConceptPlugin setParameter(Enum<PARAMETERS> key, T value) {
        return setParameter(key.toString(), value);
    }

    protected <T extends Comparable<T>> ParameterConceptPlugin setParameter(Enum<PARAMETERS> key, String value) {
        return setParameter(key.toString(), value);
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Parameter Concept";
    }

    @Override
    public String getToolTipText() {
        return "Parameter Concept Plugin";
    }

    @Override
    public String getCitation() {
        return "Casstevens T (2014) TASSEL: Self-Describing Plugins. Publication 1:1.";
    }

    private static final int TEXT_FIELD_WIDTH = 20;

    boolean parametersIsSet = true;

    /**
     * Generates dialog based on this plugins define parameters.
     *
     * @param <T>
     *
     * @return true if OK clicked, false if canceled
     */
    protected <T extends Comparable<T>> boolean setParametersViaGUI() {

        final JDialog dialog = new JDialog(getParentFrame(), null, true);

        final Map<String, JComponent> parameterFields = new HashMap<>();

        parametersIsSet = true;

        JButton okButton = new JButton();
        okButton.setActionCommand("Ok");
        okButton.setText("Ok");
        okButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    for (final PluginParameterTerry<?> current : getParameterInstances()) {
                        JComponent component = parameterFields.get(current.cmdLineName());
                        if (component instanceof JTextField) {
                            String input = ((JTextField) component).getText().trim();
                            setParameter(current.cmdLineName(), input);
                        } else if (component instanceof JCheckBox) {
                            if (((JCheckBox) component).isSelected()) {
                                setParameter(current.cmdLineName(), Boolean.TRUE);
                            } else {
                                setParameter(current.cmdLineName(), Boolean.FALSE);
                            }
                        }
                    }
                } catch (Exception ex) {
                    StringBuilder builder = new StringBuilder();
                    builder.append("Problem Setting Parameters: ");
                    builder.append("\n");
                    builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(ex), 50));
                    String str = builder.toString();
                    DialogUtils.showError(str, getParentFrame());
                    return;
                }
                dialog.setVisible(false);
            }
        });
        JButton cancelButton = new JButton();
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                parametersIsSet = false;
                dialog.setVisible(false);
            }
        });

        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));

        for (final PluginParameterTerry<?> current : getParameterInstances()) {
            if (current.valueType().isAssignableFrom(Boolean.class)) {
                JCheckBox check = new JCheckBox(current.guiName());
                if (current.value() == Boolean.TRUE) {
                    check.setSelected(true);
                } else {
                    check.setSelected(false);
                }
                JPanel temp = new JPanel(new FlowLayout(FlowLayout.LEFT));
                temp.add(check);
                panel.add(temp);
                parameterFields.put(current.cmdLineName(), check);
            } else {
                final JTextField field;
                if ((current.fileType() == PluginParameterTerry.FILE_TYPE.IN) || (current.fileType() == PluginParameterTerry.FILE_TYPE.OUT)) {
                    field = new JTextField(TEXT_FIELD_WIDTH - 5);
                } else {
                    field = new JTextField(TEXT_FIELD_WIDTH);
                }

                if (current.value() != null) {
                    field.setText(current.value().toString());
                }

                if (current.range() != null) {
                    field.addFocusListener(new FocusAdapter() {
                        @Override
                        public void focusLost(FocusEvent e) {
                            String input = field.getText().trim();
                            try {
                                T value = convert(input, (Class<T>) current.valueType());
                                if (!((Range<T>) current.range()).contains(value)) {
                                    JOptionPane.showMessageDialog(dialog, current.guiName() + " range: " + current.range().toString());
                                    field.setText(getParameterInstance(current.cmdLineName()).value().toString());
                                }
                            } catch (Exception ex) {
                                JOptionPane.showMessageDialog(dialog, current.guiName() + ": " + ex.getMessage());
                                field.setText(getParameterInstance(current.cmdLineName()).value().toString());
                            }
                        }
                    });
                }

                if (current.fileType() == PluginParameterTerry.FILE_TYPE.IN) {
                    panel.add(getLine(current.guiName(), field, getOpenFile(dialog, field)));
                } else if (current.fileType() == PluginParameterTerry.FILE_TYPE.OUT) {
                    panel.add(getLine(current.guiName(), field, getSaveFile(dialog, field)));
                } else {
                    panel.add(getLine(current.guiName(), field, null));
                }

                parameterFields.put(current.cmdLineName(), field);
            }
        }

        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.add(panel, getToolTipText());

        JPanel pnlButtons = new JPanel();
        pnlButtons.setLayout(new FlowLayout());
        pnlButtons.add(okButton);
        pnlButtons.add(cancelButton);
        dialog.getContentPane().add(tabbedPane, BorderLayout.CENTER);
        dialog.getContentPane().add(pnlButtons, BorderLayout.SOUTH);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setLocationRelativeTo(getParentFrame());
        dialog.setVisible(true);
        return parametersIsSet;

    }

    private JPanel getLine(String label, JTextField ref, JButton button) {

        JPanel result = new JPanel(new FlowLayout(FlowLayout.RIGHT));

        result.add(new JLabel(label));
        ref.setEditable(true);
        ref.setHorizontalAlignment(JTextField.LEFT);
        ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
        ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
        ref.setMaximumSize(ref.getPreferredSize());
        result.add(ref);
        if (button != null) {
            result.add(button);
        }

        return result;

    }

    private JButton getOpenFile(final JDialog parent, final JTextField textField) {

        final JFileChooser fileChooser = new JFileChooser(TasselPrefs.getOpenDir());

        JButton result = new JButton("Browse");

        result.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                if (fileChooser.showOpenDialog(parent) == JFileChooser.APPROVE_OPTION) {
                    File file = fileChooser.getSelectedFile();
                    textField.setText(file.getPath());
                    TasselPrefs.putOpenDir(fileChooser.getCurrentDirectory().getPath());
                }
            }

        });

        return result;
    }

    private JButton getSaveFile(final JDialog parent, final JTextField textField) {

        final JFileChooser fileChooser = new JFileChooser(TasselPrefs.getSaveDir());

        JButton result = new JButton("Browse");

        result.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                if (fileChooser.showSaveDialog(parent) == JFileChooser.APPROVE_OPTION) {
                    File file = fileChooser.getSelectedFile();
                    textField.setText(file.getPath());
                    TasselPrefs.putSaveDir(fileChooser.getCurrentDirectory().getPath());
                }
            }

        });

        return result;
    }

    public static class Builder {

        private final ParameterConceptPlugin myPlugin;

        public Builder(Frame frame) {
            myPlugin = new ParameterConceptPlugin(frame);
        }

        protected <T extends Comparable<T>> Builder setParameter(Enum<PARAMETERS> key, T value) {
            myPlugin.setParameter(key.toString(), value);
            return this;
        }

        public ParameterConceptPlugin build() {
            return myPlugin;
        }

    }
}
