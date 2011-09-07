/*
 * QQDisplayPlugin.java
 *
 * Created on December 13, 2010
 *
 */

package net.maizegenetics.baseplugins;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.net.URL;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.DataSet;
import java.util.List;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import net.maizegenetics.baseplugins.chart.XYScatterAndLinePanel;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.Datum;

/**
 *
 * @author yz79
 */
public class QQDisplayPlugin extends AbstractDisplayPlugin{

    /** Creates a new instance of QQDisplayPlugin */
    public QQDisplayPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {
        List <Datum> tableInList=input.getDataOfType(TableReport.class);
        if(tableInList.size()!=1) {
            String message="Invalid selection.  Please select one table result.";
            if(isInteractive()) {JOptionPane.showMessageDialog(getParentFrame(), message);} else {System.out.println(message);}
            return null;
        }
        TableReport myTableReport = (TableReport)tableInList.get(0).getData();
        if(isInteractive()) {
            try {
                QQDisplayPluginDialog myDialog = new QQDisplayPluginDialog(this, myTableReport);
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
            }
            catch(Exception ex) {
                ex.printStackTrace();
                JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create QQ plot " + ex);
            } catch(Error er) {
                er.printStackTrace();
                JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create QQ plot " + er);
            }
        }
        return null;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        return null;
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "QQ Plot";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Display QQ Plot";
    }

}

class QQDisplayPluginDialog extends JDialog {

    QQDisplayPlugin myQQDisplayPlugin;
    XYScatterAndLinePanel myQQPlot;
    TableReport myTableReport;

    JButton myCloseButton = new JButton();
    JPanel myMainPanel;
    JPanel myOptionPanel = new JPanel();
    

    public QQDisplayPluginDialog(QQDisplayPlugin plugin, TableReport theTableReport) {
        super(plugin.getParentFrame(), "QQ Plot", false);
        myQQDisplayPlugin = plugin;
        myTableReport = theTableReport;
        try {
            jbInit();
            myQQPlot = new XYScatterAndLinePanel(theTableReport);
//            myQQFigurePanel = new QQComponent(theTableReport);
            getContentPane().add(myQQPlot, BorderLayout.CENTER);
            pack();
        } catch(Exception ex) {
            ex.printStackTrace();
            JOptionPane.showMessageDialog(this.getParent(), "Unable to create QQ plot " + ex);
        }
        repaint();
    }

    void jbInit() throws Exception {
        myCloseButton.setText("Close");
        myCloseButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                closeButton_actionPerformed(e);
            }
        });

        myOptionPanel.add(myCloseButton, new GridBagConstraints(0, 0, 0, 0, 0.0, 0.0
                ,GridBagConstraints.SOUTHEAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        getContentPane().add(myOptionPanel);
    }

    void closeButton_actionPerformed(ActionEvent e) {
        dispose();
    }

}