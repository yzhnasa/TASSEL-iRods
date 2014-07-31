/*
 * AbstractPlugin.java
 *
 * Created on December 22, 2006, 5:03 PM
 *
 */
package net.maizegenetics.plugindef;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

import java.io.File;
import java.lang.reflect.Field;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
abstract public class AbstractPlugin implements Plugin {

    private static final Logger myLogger = Logger.getLogger(AbstractPlugin.class);

    private final List<PluginListener> myListeners = new ArrayList<>();
    private final List<Plugin> myInputs = new ArrayList<>();
    private final Frame myParentFrame;
    private final boolean myIsInteractive;
    private boolean myTrace = false;
    private boolean myThreaded = false;

    /**
     * Creates a new instance of AbstractPlugin
     */
    public AbstractPlugin() {
        this(null, true);
    }

    /**
     * Creates a new instance of AbstractPlugin
     */
    public AbstractPlugin(Frame parentFrame, boolean isInteractive) {
        myParentFrame = parentFrame;
        myIsInteractive = isInteractive;
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

    @Override
    public DataSet processData(DataSet input) {
        throw new UnsupportedOperationException();
    }

    private List<PluginParameter<?>> getParameterInstances() {

        List<PluginParameter<?>> result = new ArrayList<>();
        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameter.class)) {
                current.setAccessible(true);
                try {
                    PluginParameter<?> parameter = (PluginParameter) current.get(this);
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
            if (current.getType().isAssignableFrom(PluginParameter.class)) {
                try {
                    current.setAccessible(true);
                    PluginParameter<?> parameter = (PluginParameter) current.get(this);
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

    private PluginParameter<?> getParameterInstance(String key) {
        try {
            Field field = getParameterField(key);
            return (PluginParameter) field.get(this);
        } catch (Exception e) {
            return null;
        }
    }

    public static <T> T convert(String input, Class<T> outputClass) {
        try {
            if (outputClass.isEnum()) {
                return (T) Enum.valueOf(outputClass.asSubclass(Enum.class), input);
            } else {
                return input == null ? null : outputClass.getConstructor(String.class).newInstance(input);
            }
        } catch (Exception nfe) {
            throw new IllegalArgumentException("convert: Problem converting: " + input + " to " + outputClass.getName());
        }
    }

    @Override
    public void setParameters(String[] args) {

        if (args != null) {

            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (args[i].startsWith("-")) {
                    arg = arg.substring(1);
                    PluginParameter<?> parameter = getParameterInstance(arg);
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

        for (PluginParameter<?> current : getParameterInstances()) {

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

            if (current.fileType() == PluginParameter.FILE_TYPE.IN) {
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

            if (current.fileType() == PluginParameter.FILE_TYPE.OUT) {
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

    protected void printParameterValues() {
        List<PluginParameter<?>> parameters = getParameterInstances();
        if ((parameters == null) || (parameters.size() == 0)) {
            return;
        }
        StringBuilder builder = new StringBuilder();
        builder.append("\n");
        builder.append(Utils.getBasename(getClass().getName()));
        builder.append(" Parameters\n");
        for (PluginParameter<?> current : parameters) {
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
        for (PluginParameter<?> current : getParameterInstances()) {
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
                if (current.valueType().isEnum()) {
                    builder.append(" [");
                    Comparable[] values = current.valueType().getEnumConstants();
                    for (int i = 0; i < values.length; i++) {
                        if (i != 0) {
                            builder.append(" ");
                        }
                        builder.append(values[i].toString());
                    }
                    builder.append("]");
                } else {
                    builder.append(" ");
                    builder.append(current.range().toString());
                }
            }
            if (current.defaultValue() != null) {
                builder.append(" (Default: ");
                builder.append(current.defaultValue());
                builder.append(")");
            }
            if (current.required()) {
                builder.append(" (required)");
            }
            builder.append("\n");
        }
        myLogger.info(builder.toString());
    }

    @Override
    public Comparable getParameter(Enum key) {
        return getParameterInstance(key.toString()).value();
    }

    @Override
    public Comparable getParameter(String key) {
        return getParameterInstance(key).value();
    }

    @Override
    public Plugin setParameter(String key, Comparable value) {

        PluginParameter parameter = null;
        try {

            Field field = getParameterField(key);
            parameter = (PluginParameter) field.get(this);
            if (parameter == null) {
                throw new IllegalArgumentException("setParameter: Unknown Parameter: " + key);
            }
            if ((parameter.range() != null) && (!parameter.range().contains(value))) {
                throw new IllegalArgumentException("setParameter: " + parameter.cmdLineName() + " value: " + value.toString() + " outside range: " + parameter.range().toString());
            }
            PluginParameter newParameter = new PluginParameter<>(parameter, value);
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

    @Override
    public Plugin setParameter(String key, String value) {

        PluginParameter parameter = null;
        try {
            parameter = getParameterInstance(key);
            return setParameter(key, (Comparable) convert(value, parameter.valueType()));
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

    @Override
    final public Plugin setParameter(Enum key, Comparable value) {
        return setParameter(key.toString(), value);
    }

    @Override
    final public Plugin setParameter(Enum key, String value) {
        return setParameter(key.toString(), value);
    }

    private static final int TEXT_FIELD_WIDTH = 25;

    private boolean parametersAreSet = true;

    /**
     * Generates dialog based on this plugins define parameters.
     *
     * @param <T>
     *
     * @return true if OK clicked, false if canceled
     */
    private <T extends Comparable<T>> boolean setParametersViaGUI() {

        final JDialog dialog = new JDialog(getParentFrame(), null, true);

        final Map<String, JComponent> parameterFields = new HashMap<>();

        parametersAreSet = true;

        JButton okButton = new JButton();
        okButton.setActionCommand("Ok");
        okButton.setText("Ok");
        okButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    for (final PluginParameter<?> current : getParameterInstances()) {
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
                        } else if (component instanceof JComboBox) {
                            Enum temp = (Enum) ((JComboBox) component).getSelectedItem();
                            setParameter(current.cmdLineName(), temp);
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
                parametersAreSet = false;
                dialog.setVisible(false);
            }
        });

        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));

        for (final PluginParameter<?> current : getParameterInstances()) {
            if (current.valueType().isEnum()) {
                JComboBox menu = new JComboBox();
                Comparable[] values = current.valueType().getEnumConstants();
                for (Comparable item : values) {
                    menu.addItem(item);
                }
                menu.setSelectedItem(current.value());
                JPanel temp = new JPanel(new FlowLayout(FlowLayout.CENTER));
                temp.add(new JLabel(current.guiName()));
                temp.add(menu);
                panel.add(temp);
                parameterFields.put(current.cmdLineName(), menu);
            } else if (current.valueType().isAssignableFrom(Boolean.class)) {
                JCheckBox check = new JCheckBox(current.guiName());
                if (current.value() == Boolean.TRUE) {
                    check.setSelected(true);
                } else {
                    check.setSelected(false);
                }
                JPanel temp = new JPanel(new FlowLayout(FlowLayout.CENTER));
                temp.add(check);
                panel.add(temp);
                parameterFields.put(current.cmdLineName(), check);
            } else {
                final JTextField field;
                if ((current.fileType() == PluginParameter.FILE_TYPE.IN) || (current.fileType() == PluginParameter.FILE_TYPE.OUT)) {
                    field = new JTextField(TEXT_FIELD_WIDTH - 8);
                } else {
                    field = new JTextField(TEXT_FIELD_WIDTH);
                }

                if (current.value() != null) {
                    field.setText(current.value().toString());
                }

                field.addFocusListener(new FocusAdapter() {
                    @Override
                    public void focusLost(FocusEvent e) {
                        String input = field.getText().trim();
                        try {
                            if (!current.acceptsValue(input)) {
                                JOptionPane.showMessageDialog(dialog, current.guiName() + " range: " + current.range().toString());
                                field.setText(getParameterInstance(current.cmdLineName()).value().toString());
                            }
                        } catch (Exception ex) {
                            JOptionPane.showMessageDialog(dialog, current.guiName() + ": " + ex.getMessage());
                            field.setText(getParameterInstance(current.cmdLineName()).value().toString());
                        }
                    }
                });

                if (current.fileType() == PluginParameter.FILE_TYPE.IN) {
                    panel.add(getLine(current.guiName(), field, getOpenFile(dialog, field)));
                } else if (current.fileType() == PluginParameter.FILE_TYPE.OUT) {
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
        return parametersAreSet;

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

    /**
     * Returns menu that can be added to main menu bar.
     *
     * @return menu
     */
    @Override
    public JMenu getMenu() {
        return null;
    }

    /**
     * Sets up this plugin to receive input from another plugin.
     *
     * @param input input
     */
    @Override
    public void receiveInput(Plugin input) {

        if (input == null) {
            throw new IllegalArgumentException("AbstractPlugin: receiveInput: input can not be null.");
        }

        if (!myInputs.contains(input)) {
            myInputs.add(input);
        }

        input.addListener(this);

    }

    /**
     * GUI Panel for this plugin.
     *
     * @return panel
     */
    @Override
    public JPanel getPanel() {
        return null;
    }

    /**
     * If interactive = true, the plugin will create dialogs and panels to
     * interacts with the user
     *
     * @return boolean
     */
    @Override
    public boolean isInteractive() {
        return myIsInteractive;
    }

    /**
     * Parent Frame for this plugin. Can be null.
     *
     * @return frame
     */
    @Override
    public Frame getParentFrame() {
        return myParentFrame;
    }

    /**
     * Adds listener to this plugin.
     *
     * @param listener listener to add
     */
    @Override
    public void addListener(PluginListener listener) {

        synchronized (myListeners) {
            if ((listener != null) && (!myListeners.contains(listener))) {
                myListeners.add(listener);
            }
        }

    }

    protected List<PluginListener> getListeners() {
        return myListeners;
    }

    public List<Plugin> getInputs() {
        return myInputs;
    }

    /**
     * Returns data set after complete.
     *
     * @param event event
     */
    public void fireDataSetReturned(PluginEvent event) {

        synchronized (myListeners) {
            Iterator<PluginListener> itr = myListeners.iterator();
            while (itr.hasNext()) {
                try {
                    if (myThreaded) {
                        PluginListener current = itr.next();
                        ThreadedPluginListener thread = new ThreadedPluginListener(current, event);
                        thread.start();
                    } else {
                        PluginListener current = itr.next();
                        current.dataSetReturned(event);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

    }

    /**
     * Returns data set after complete.
     *
     * @param data data set
     */
    protected void fireDataSetReturned(DataSet data) {
        fireDataSetReturned(new PluginEvent(data));
    }

    /**
     * Returns progress of execution.
     *
     * @param event event
     */
    protected void fireProgress(PluginEvent event) {

        synchronized (myListeners) {
            Iterator<PluginListener> itr = myListeners.iterator();
            while (itr.hasNext()) {
                PluginListener current = itr.next();
                current.progress(event);
            }
        }

        DataSet ds = (DataSet) event.getSource();
        if (ds != null) {
            List percentage = ds.getDataOfType(Integer.class);

            if (percentage.size() > 0) {
                Datum datum = (Datum) percentage.get(0);
                Integer percent = (Integer) datum.getData();
                if (percent == 100) {
                    myLogger.info(getClass().getName() + "  Citation: " + getCitation());
                }
            }
        }

    }

    /**
     * Returns progress of execution.
     *
     * @param percent percentage between 0 and 100 inclusive.
     */
    protected void fireProgress(Integer percent) {

        if ((percent < 0) || (percent > 100)) {
            throw new IllegalArgumentException("AbstractPlugin: fireProgress: percent must be between 0 and 100 inclusive.  arg: " + percent);
        }

        Datum percentage = new Datum("Percent", percent, null);
        fireProgress(new PluginEvent(new DataSet(percentage, this)));

    }

    @Override
    public String getCitation() {
        return "Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. (2007) TASSEL: Software for association mapping of complex traits in diverse samples. Bioinformatics 23:2633-2635.";
    }

    //
    // Methods for PluginListener.
    //
    /**
     * Returns data set after complete.
     *
     * @param event event
     */
    @Override
    public void dataSetReturned(PluginEvent event) {

        DataSet input = (DataSet) event.getSource();

        performFunction(input);

    }

    /**
     * No operation for this abstract class.
     */
    @Override
    public void progress(PluginEvent event) {
        // The default action of a plugin is to do
        // nothing when another plugin reports its
        // progress.   This is intended to be implemented
        // by GUI applications to show the user the
        // progress of an interactive action.
    }

    public void reverseTrace(int indent) {

        if (myTrace) {
            return;
        }

        indent(indent);
        System.out.println(getClass().getName());

        Iterator<Plugin> itr = myInputs.iterator();
        while (itr.hasNext()) {
            try {
                AbstractPlugin current = (AbstractPlugin) itr.next();
                current.reverseTrace(indent + 3);
            } catch (Exception e) {
                // do nothing
            }
        }

        myTrace = true;

    }

    public void trace(int indent) {

        if (myTrace) {
            return;
        }

        indent(indent);
        System.out.println(getClass().getName());

        Iterator<PluginListener> itr = myListeners.iterator();
        while (itr.hasNext()) {
            try {
                AbstractPlugin current = (AbstractPlugin) itr.next();
                current.trace(indent + 3);
            } catch (Exception e) {
                // do nothing
            }
        }

        myTrace = true;

    }

    private void indent(int indent) {

        for (int i = 0; i < indent; i++) {
            System.out.print(" ");
        }

    }

    @Override
    public void setThreaded(boolean threaded) {
        myThreaded = threaded;
    }

    @Override
    public boolean cancel() {
        return false;
    }

    @Override
    public void run() {
        performFunction(null);
    }

    @Override
    public void progress(int percent, Object meta) {
        fireProgress(percent);
    }

}
