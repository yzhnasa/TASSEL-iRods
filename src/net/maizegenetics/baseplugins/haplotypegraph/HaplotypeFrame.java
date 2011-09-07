package net.maizegenetics.baseplugins.haplotypegraph;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.tassel.DataTreePanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jul 6, 2004
 * Time: 8:01:19 AM
 * To change this template use File | Settings | File Templates.
 */
public class HaplotypeFrame extends JFrame {

    HaplotypePanel haplotypePanel;
    JPanel southPanel;
    JRadioButton minSpanTree;
    JRadioButton minSpanTreeSame;
    JRadioButton minSpanTreeThreshold;
    JButton svgSave;
    JLabel thresholdLabel;
    JTextField thresholdValue;
    ButtonGroup buttonGroup;
    DataTreePanel theDataTreePanel;
    JungTestPanel jtp;

    // this constructor is for testing, in practice and once integrated it will most likely take in an AnnotationAlignment
    public HaplotypeFrame(File file) {
        super();
        haplotypePanel = new HaplotypePanel(file);
        jtp = new JungTestPanel();
        init();
    }

    public HaplotypeFrame(DataTreePanel theDataTreePanel, Alignment aa) {
        super();
        haplotypePanel = new HaplotypePanel(aa);
        jtp = new JungTestPanel(aa, this);
        this.theDataTreePanel = theDataTreePanel;
        init();
    }

    public void run() {
        while (true) {
            repaint();
        }
    }

    public DataTreePanel getDataTreePanel() {
        return theDataTreePanel;
    }

    private void init() {
        // initialize everything
        southPanel = new JPanel();
        minSpanTree = new JRadioButton();
        minSpanTreeSame = new JRadioButton();
        minSpanTreeThreshold = new JRadioButton();
        svgSave = new JButton();
        thresholdLabel = new JLabel();
        thresholdValue = new JTextField();
        buttonGroup = new ButtonGroup();

        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();

        // setup the buttons
        setUpButtons();

        // set up the haplotypePanel
        haplotypePanel.setBackground(Color.white);

        // set up the theshold fields
        thresholdLabel.setText("Enter the threshold:");
        thresholdValue.setText("0");
        thresholdValue.setColumns(10);
        thresholdValue.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                thresholdValue_actionPerformed(e);
            }
        });

        // set up the southPanel
        southPanel.add(minSpanTree);
        southPanel.add(minSpanTreeSame);
        southPanel.add(minSpanTreeThreshold);
        southPanel.add(thresholdLabel);
        southPanel.add(thresholdValue);
        southPanel.add(svgSave);



        // set up the frame
        this.setTitle("Haplotype Viewer");

        this.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);

        this.setSize(screenSize.width * 3 / 4 + 30, screenSize.height * 3 / 4 + 30);

        Container c = this.getContentPane();
        Box mainBox = new Box(BoxLayout.Y_AXIS);
        c.add(mainBox);
        mainBox.add(jtp);

        this.show();
    }

    private void setUpButtons() {
        minSpanTree.setText("Minimum Connection");
        minSpanTree.setToolTipText("view the minimum spanning tree connection of the haplotyes");
        minSpanTree.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                minSpanTree_actionPerformed(e);
            }
        });

        minSpanTreeSame.setText("Default View");
        minSpanTreeSame.setToolTipText("view the minimum spanning tree connection where all edges connecting a haplotype to the tree with he same difference are shown");
        minSpanTreeSame.setSelected(true);
        minSpanTreeSame.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                minSpanTreeSame_actionPerformed(e);
            }
        });

        minSpanTreeThreshold.setText("View Threshold Edges");
        minSpanTreeThreshold.setToolTipText("view the connections where the difference between haplotypes is withing the given threshold");
        minSpanTreeThreshold.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                minSpanTreeThreshold_actionPerformed(e);
            }
        });

        svgSave.setText("Save");
        svgSave.setToolTipText("save as a scalable vector graphic");
        svgSave.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                svgSave_actionPerformed(e);
            }
        });

        buttonGroup.add(minSpanTree);
        buttonGroup.add(minSpanTreeSame);
        buttonGroup.add(minSpanTreeThreshold);
    }

    private void minSpanTree_actionPerformed(ActionEvent e) {
        haplotypePanel.setDrawType(HaplotypePanel.MINIMUM_SPANNING_TREE);
        repaint();
    }

    private void minSpanTreeSame_actionPerformed(ActionEvent e) {
        haplotypePanel.setDrawType(HaplotypePanel.MINIMUM_SPANNING_TREE_SAME);
        repaint();
    }

    private void minSpanTreeThreshold_actionPerformed(ActionEvent e) {
        int threshold = 0;
        threshold = Integer.parseInt(thresholdValue.getText());
        haplotypePanel.setDrawType(HaplotypePanel.MINIMUM_SPANNING_TREE_THRESHOLD);
        haplotypePanel.setDrawThreshold(threshold);
        repaint();
    }

    private void thresholdValue_actionPerformed(ActionEvent e) {
        if (minSpanTreeThreshold.isSelected()) {
            minSpanTreeThreshold_actionPerformed(e);
            thresholdValue.selectAll();
        }
    }

    private void svgSave_actionPerformed(ActionEvent e) {
        JFileChooser fc = new JFileChooser();
        int returnVal = fc.showSaveDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION);
        {
            File file = fc.getSelectedFile();
            haplotypePanel.svgSave(file);
        }
    }
}
