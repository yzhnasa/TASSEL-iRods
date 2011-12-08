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
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.DataSet;
import java.util.List;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
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
                PlotOptionsQQDialog myOptions = new PlotOptionsQQDialog(this.getParentFrame(), getTraits(myTableReport), splitTable(myTableReport));
                myOptions.setLocationRelativeTo(getParentFrame());
                myOptions.setVisible(true);
                if (myOptions.isCanceled() == false) {
                    QQDisplayPluginDialog myDialog = new QQDisplayPluginDialog(this.getParentFrame(), myTableReport, myOptions.getSliderValue());
                    myDialog.setLocationRelativeTo(getParentFrame());
                    myDialog.setVisible(true);
                }
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

    private ArrayList<Integer> splitTable(TableReport table) {
        ArrayList<Integer> indexes = new ArrayList<Integer>();
        int numRows = table.getRowCount();
        String previousTrait = "";
        for (int i = 0; i < numRows; i++) {
            if (!previousTrait.equals((String)table.getValueAt(i, 0))) {
                if (!((String)table.getValueAt(i, 1)).equals("None")) {
                    indexes.add(new Integer(i));
                    previousTrait = (String)table.getValueAt(i, 0);
                    if (i > 1) {
                        indexes.add(new Integer(i));
                    }
                } else if (i != 0) {
                    indexes.add(new Integer(i));
                    indexes.add(new Integer(i+1));
                    previousTrait = (String)table.getValueAt(i + 1, 0);
                }
            }
        }
        indexes.add(new Integer(numRows));
        return indexes;
    }

    private String[] getTraits(TableReport table) {
        ArrayList<String> traitArray = new ArrayList<String>();
        int numRows = table.getRowCount();
        String previousTrait = "";
        for (int i = 0; i < numRows; i++) {
            if (!previousTrait.equals((String)table.getValueAt(i, 0))) {
                previousTrait = (String)table.getValueAt(i, 0);
                traitArray.add(previousTrait);
            }
        }

        String[] traits = new String[traitArray.size()];
        for (int i = 0; i < traitArray.size(); i++) {
            traits[i] = traitArray.get(i);
        }
        return traits;
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

    XYScatterAndLinePanel myQQPlot;
    TableReport myTableReport;

    JButton myCloseButton = new JButton();
    JPanel myMainPanel;
    JPanel myOptionPanel = new JPanel();
    

    public QQDisplayPluginDialog(Frame f, TableReport theTableReport, int countToDisplay) {
        super(f, "QQ Plot", false);
        myTableReport = theTableReport;
        try {
            jbInit();
            myQQPlot = new XYScatterAndLinePanel(theTableReport, countToDisplay);
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

class PlotOptionsQQDialog extends JDialog {

    boolean isCanceled = true;
    private JPanel mainPanel = new JPanel();
    private JButton okayButton = new JButton();
    private JButton cancelButton = new JButton();
    private JSlider slider = new JSlider();
    private JLabel countLabel1 = new JLabel();
    private JLabel countLabel2 = new JLabel();
    private JTextField countTextField = new JTextField();
    private GridBagLayout gridBagLayout2 = new GridBagLayout();

    public PlotOptionsQQDialog(Frame f, String[] traits, ArrayList<Integer> indexes) {
        super(f, "QQ Plot Options", true);

        int numSites = indexes.get(1) - indexes.get(0);
        slider.setMinimum(1);
        slider.setMaximum(numSites);
        slider.setValue((int)(numSites * 0.05));
        countTextField.setText("" + (int)(numSites * 0.05));

        try {
            initUI();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void initUI() throws Exception {

        // mainPanel.setBackground(SystemColor.menu);
        mainPanel.setMinimumSize(new Dimension(400,80));
        mainPanel.setPreferredSize(new Dimension(400, 80));
        mainPanel.setLayout(gridBagLayout2);

        GridBagConstraints c = new GridBagConstraints();

        slider.addChangeListener(new javax.swing.event.ChangeListener() {

            public void stateChanged(ChangeEvent ce) {
                slider_actionPerformed(ce);
            }
        });
        c.gridx = 0;
        c.gridy = 0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridwidth = 3;
        mainPanel.add(slider, c);

        countLabel1.setText("Plotting");
        c.gridx = 0;
        c.gridy = 1;
        c.gridwidth = 1;
        mainPanel.add(countLabel1, c);

        countTextField.addKeyListener(new java.awt.event.KeyListener() {

            public void keyTyped(KeyEvent ke) {
//                countTextField_keyTyped(ke);
            }

            public void keyPressed(KeyEvent ke) {
//                throw new UnsupportedOperationException("Not supported yet.");
            }

            public void keyReleased(KeyEvent ke) {
                countTextField_keyTyped(ke);
            }
        });

        c.gridx = 1;
        c.gridy = 1;
        c.weighty = 1;
        mainPanel.add(countTextField, c);

        countLabel2.setText("most significant out of " + slider.getMaximum() + " per trait.");
        c.gridx = 2;
        c.gridy = 1;
        c.weighty = 0;
        mainPanel.add(countLabel2, c);

        okayButton.setMaximumSize(new Dimension(63, 27));
        okayButton.setMinimumSize(new Dimension(63, 27));
        okayButton.setText("Okay");
        okayButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                okayButton_actionPerformed(e);
            }
        });
        c.gridx = 1;
        c.gridy = 2;
        mainPanel.add(okayButton, c);

        cancelButton.setMaximumSize(new Dimension(63, 27));
        cancelButton.setMinimumSize(new Dimension(63, 27));
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });
        c.gridx = 2;
        c.gridy = 2;
        mainPanel.add(cancelButton, c);

        this.add(mainPanel, BorderLayout.CENTER);
    }

    private void okayButton_actionPerformed(ActionEvent e) {
        isCanceled = false;
        setVisible(false);
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        isCanceled = true;
        setVisible(false);
    }

    private void slider_actionPerformed(ChangeEvent cd) {
        countTextField.setText("" + slider.getValue());
    }

    private void countTextField_keyTyped(KeyEvent e) {
        try {
            if (!countTextField.getText().equals("")) {
                int value = Integer.valueOf(countTextField.getText());
                if (value >= slider.getMinimum() && value <= slider.getMaximum()) {
                    slider.setValue(value);
                } else if (value <= slider.getMinimum()) {
                    slider.setValue(slider.getMinimum());
                    countTextField.setText("" + slider.getMinimum());
                } else if (value >= slider.getMaximum()) {
                    slider.setValue(slider.getMaximum());
                    countTextField.setText("" + slider.getMaximum());
                }
            }
        } catch (NumberFormatException nfe) {
            countTextField.setText("" + slider.getValue());
        }
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    public int getSliderValue() {
        return slider.getValue();
    }

    public Dimension getMinimumSize() {
        return new Dimension(600, 600);
    }
}