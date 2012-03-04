/*
 * ConvertSBitTBitPlugin
 */
package net.maizegenetics.baseplugins;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTabbedPane;

import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.SBitAlignment;
import net.maizegenetics.pal.alignment.TBitAlignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class ConvertSBitTBitPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ConvertSBitTBitPlugin.class);

    public enum CONVERT_TYPE {

        sbit, tbit
    };
    public CONVERT_TYPE myType = CONVERT_TYPE.sbit;

    public ConvertSBitTBitPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        String name = null;

        try {

            List inputData = input.getDataOfType(Alignment.class);
            if (inputData.size() != 1) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Please Select a Single Alignment.");
                } else {
                    myLogger.error("performFunction: Please Select a Single Alignment.");
                }
                return null;
            }

            Datum alignDatum = (Datum) inputData.get(0);
            Alignment alignment = (Alignment) alignDatum.getData();
            name = alignDatum.getName();

            if (isInteractive()) {
                ConvertSBitTBitPluginDialog theDialog = new ConvertSBitTBitPluginDialog();
                theDialog.setLocationRelativeTo(getParentFrame());
                theDialog.setVisible(true);
                if (theDialog.isCancel()) {
                    return null;
                }
                myType = theDialog.getType();
                theDialog.dispose();
            }

            Alignment newAlignment = convertAlignment(alignment);
            if (newAlignment == null) {
                showErrorMessage(name, null);
                return null;
            } else if (alignment == newAlignment) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Nothing To Change");
                    return null;
                }
            }
            DataSet tds = new DataSet(new Datum(alignDatum.getName() + "_" + myType, newAlignment, null), this);
            fireDataSetReturned(new PluginEvent(tds, ConvertSBitTBitPlugin.class));
            return tds;

        } catch (Exception e) {
            e.printStackTrace();
            showErrorMessage(name, e);
            return null;
        } finally {
            fireProgress(100);
        }
    }

    private Alignment convertAlignment(Alignment alignment) {
        Alignment result = null;
        if (myType == CONVERT_TYPE.sbit) {
            return SBitAlignment.getInstance(alignment);
        } else if (myType == CONVERT_TYPE.tbit) {
            return TBitAlignment.getInstance(alignment);
        } else {
            throw new IllegalStateException("ConvertSBitTBitPlugin: convertAlignment: Unknown type: " + myType);
        }
    }

    private void showErrorMessage(String name, Exception e) {
        StringBuilder builder = new StringBuilder();
        builder.append("Error converting: ");
        builder.append(name);
        builder.append(" to: ");
        builder.append(myType);
        if (e != null) {
            builder.append("\n");
            builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
        }
        String str = builder.toString();

        if (isInteractive()) {
            DialogUtils.showError(str, getParentFrame());
        } else {
            myLogger.error(str);
        }
    }

    public ImageIcon getIcon() {
        return null;
    }

    public String getButtonName() {
        return "TtoS";
    }

    public String getToolTipText() {
        return "Converts Between Taxa and Site Optimized Data Sets.";
    }

    public CONVERT_TYPE getType() {
        return myType;
    }

    public void setType(CONVERT_TYPE type) {
        myType = type;
    }
}

class ConvertSBitTBitPluginDialog extends JDialog {

    private ConvertSBitTBitPlugin.CONVERT_TYPE myType = ConvertSBitTBitPlugin.CONVERT_TYPE.sbit;
    private JTabbedPane myTabbedPane = new JTabbedPane();
    private ButtonGroup myButtonGroup = new ButtonGroup();
    private JRadioButton mySBit = new JRadioButton();
    private JRadioButton myTBit = new JRadioButton();
    private boolean myIsCancel = true;

    public ConvertSBitTBitPluginDialog() {
        super((Frame) null, "Converts Between Taxa and Site Optimized Data Sets.", true);


        JButton okButton = new JButton();
        okButton.setActionCommand("Ok");
        okButton.setText("Ok");
        okButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                myIsCancel = false;
                setVisible(false);
            }
        });
        JButton closeButton = new JButton();
        closeButton.setText("Close");
        closeButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                myIsCancel = true;
                setVisible(false);
            }
        });

        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));

        //Radio Buttons
        mySBit.setText("Optimize For Site Operations");
        myTBit.setText("Optimize For Taxa Operations");

        myButtonGroup.add(mySBit);
        myButtonGroup.add(myTBit);
        myButtonGroup.setSelected(mySBit.getModel(), true);

        panel.add(mySBit);
        panel.add(myTBit);

        myTabbedPane.add(panel, "Optimized Data Set");

        JPanel pnlButtons = new JPanel();
        pnlButtons.setLayout(new FlowLayout());
        pnlButtons.add(okButton);
        pnlButtons.add(closeButton);
        getContentPane().add(myTabbedPane, BorderLayout.CENTER);
        getContentPane().add(pnlButtons, BorderLayout.SOUTH);

        pack();
    }

    public boolean isCancel() {
        return myIsCancel;
    }

    public ConvertSBitTBitPlugin.CONVERT_TYPE getType() {
        if (mySBit.isSelected()) {
            return ConvertSBitTBitPlugin.CONVERT_TYPE.sbit;
        } else if (myTBit.isSelected()) {
            return ConvertSBitTBitPlugin.CONVERT_TYPE.tbit;
        } else {
            throw new IllegalStateException("ConvertSBitTBitPluginDialog: getType: Unknown state.");
        }
    }
}
