package net.maizegenetics.baseplugins;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.stats.MLM.Kinship;

import javax.swing.*;
import java.net.URL;
import java.awt.Container;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Author: Zhiwu Zhang
 * Date: Apr 29, 2007
 */
public class KinshipPlugin extends AbstractPlugin {

    private boolean recognizeHets = true;
    private boolean rescaleKinship = true;

    public KinshipPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);

    }

    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataSet();

            if (alignInList.isEmpty()) {
                String message = "Nothing selected. Please select pedigree data.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }

            List result = new ArrayList();
            Iterator itr = alignInList.iterator();
            while (itr.hasNext()) {

                Datum current = (Datum) itr.next();
                String datasetName = current.getName();
                Kinship kin = null;

                try {

                    if (current.getData() instanceof Alignment) {
                        //this section implements additional options for calculating kinship
                        if (isInteractive()) {
                            KinshipDialog kd = new KinshipDialog(getParentFrame());
                            kd.setLocationRelativeTo(getParentFrame());
                            kd.setVisible(true);
                            if (!kd.run) {
                                return null;

                            }
                            recognizeHets = kd.btnRelated.isSelected();
                            rescaleKinship = kd.chkRescale.isSelected();
                            kd.dispose();
                            getParentFrame().repaint();
                        }
                        Alignment theAlignment = (Alignment) current.getData();
                        kin = new Kinship(theAlignment, recognizeHets, rescaleKinship);

                    } else if (current.getData() instanceof SimplePhenotype) { //pedigree data
                        SimplePhenotype ped = (SimplePhenotype) current.getData();
                        kin = new Kinship(ped);
                    } else {
                        String message = "Invalid selection. Can't create kinship matrix from: " + datasetName;
                        if (isInteractive()) {
                            JOptionPane.showMessageDialog(getParentFrame(), message);
                        } else {
                            System.out.println(message);
                        }
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    String message = "Problem creating kinship matrix from: " + datasetName + "\n" + e.getClass().getName() + ": " + e.getMessage();
                    if (isInteractive()) {
                        JOptionPane.showMessageDialog(getParentFrame(), message);
                    } else {
                        System.out.println(message);
                        e.printStackTrace();
                    }
                }

                if (kin != null) {
                    //add kin to datatree;
                    DataSet ds = new DataSet(new Datum("kin_" + datasetName, kin.getDm(), "kinship matrix created from " + datasetName), this);
                    result.add(ds);
                    fireDataSetReturned(new PluginEvent(ds, KinshipPlugin.class));
                }

            }

            return DataSet.getDataSet(result, this);

        } finally {
            fireProgress(100);
        }

    }

    public ImageIcon getIcon() {
        URL imageURL = KinshipPlugin.class.getResource("images/Kin.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "Kinship";
    }

    public String getToolTipText() {
        return "Calculate kinship from marker data";
    }

    public void setRecognizeHets(boolean recognizeHets) {
        this.recognizeHets = recognizeHets;
    }

    public void setRescaleKinship(boolean rescaleKinship) {
        this.rescaleKinship = rescaleKinship;
    }

    class KinshipDialog extends JDialog {

        JRadioButton btnRelated = new JRadioButton("Related to homozygotes", true);
        JRadioButton btnUnrelated = new JRadioButton("Independent allele state", false);
        JCheckBox chkRescale = new JCheckBox("Rescale results between 2 and 0", false);
        boolean run = false;

        KinshipDialog(Frame parent) {
            super(parent, true);
            setLocationRelativeTo(parent);
            setTitle("Kinship Options");
            Container contents = getContentPane();

            contents.setLayout(new GridBagLayout());
            GridBagConstraints gbc = new GridBagConstraints();

            JPanel hetPanel = new JPanel();
            hetPanel.setLayout(new BoxLayout(hetPanel, BoxLayout.Y_AXIS));
            hetPanel.setBorder(BorderFactory.createTitledBorder("Model heterozygotes as"));
            ButtonGroup bg = new ButtonGroup();
            bg.add(btnRelated);
            bg.add(btnUnrelated);
            hetPanel.add(btnRelated);
            hetPanel.add(btnUnrelated);

            JButton btnCancel = new JButton("Cancel");
            btnCancel.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent arg0) {
                    run = false;
                    setVisible(false);
                }
            });

            JButton btnRun = new JButton("Run");
            btnRun.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent arg0) {
                    run = true;
                    setVisible(false);
                }
            });

            JButton btnHelp = new JButton("Help");
            btnHelp.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent arg0) {
                    String msg = "Modeling heterozygotes as related to homozygotes calculates P(IBS)= 0.5 between \n"
                            + "heterozygous loci and between a heterozygote and a homozygote. \n"
                            + "Modeling heterozygotes as an independent state uses the older version \n"
                            + "of the TASSEL kinship function, which calculates P(IBS) = 1 for a heterozygote with itself \n"
                            + "and 0 with the homozygotes.\n\n"
                            + "Choosing to model heterozygotes as an independent state and rescaling results\n"
                            + " in the same kinship matrix produced by TASSEL before the new options were added. ";
                    String myTitle = "Help for Kinship Options";
                    JOptionPane.showMessageDialog(null, msg, myTitle, JOptionPane.INFORMATION_MESSAGE);
                }
            });

            JPanel buttonBox = new JPanel();
            buttonBox.setLayout(new BoxLayout(buttonBox, BoxLayout.X_AXIS));
            buttonBox.add(btnRun);
            buttonBox.add(btnCancel);
            buttonBox.add(btnHelp);

            gbc.gridx = 0;
            gbc.gridy = 0;
            gbc.insets = new Insets(20, 10, 5, 10);
            contents.add(hetPanel, gbc);

            gbc.gridy++;
            gbc.insets = new Insets(5, 10, 5, 10);
            contents.add(chkRescale, gbc);

            gbc.gridy++;
            gbc.insets = new Insets(5, 10, 20, 10);
            contents.add(buttonBox, gbc);

            pack();
        }
    }
}
