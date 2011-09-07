package net.maizegenetics.baseplugins.haplotypegraph;

import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.GLM.ReportWriter;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Enumeration;
import java.util.LinkedList;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Aug 19, 2004
 * Time: 4:44:56 PM
 */
public class phenotypeChooserDialog extends JDialog {
    //fields
    private JList               reportList;
    private Datum            selectedBook;
    private JButton             OKButton;
    private JButton             cancelButton;
    private JTextField          textSelectedReport;
    private boolean             okButtonPressed = false;
    private HaplotypeFrame               theHapFrame;

    //constructors
    public phenotypeChooserDialog(HaplotypeFrame theHapFrame) {

        super(theHapFrame, "Choose Haplotype Effects", true);
        this.theHapFrame = theHapFrame;
        init();
    }

    //methods
    private void init() {
        this.setSize(300, 400);
        Container contentPane = getContentPane();
        contentPane.setLayout(new GridBagLayout());

        //get the haplotype reports hanging on the data tree and add them to an Object array
        Enumeration enumReports = getNode("Association").children();
        LinkedList reports = new LinkedList();
        while (enumReports.hasMoreElements()) {
            DefaultMutableTreeNode nodeReport = (DefaultMutableTreeNode) enumReports.nextElement();
            reports.add(nodeReport.getUserObject());
        }

        //create a Jlist to hold the values
        JScrollPane scrollPane;

        if (reports.size() > 0) {
            reportList = new JList(reports.toArray());
            reportList.addListSelectionListener(new ListSelectionListener() {
                public void valueChanged(ListSelectionEvent e) {
                    reportList_ValueChanged();
                }
            });
            scrollPane = new JScrollPane(reportList);
        }
        else {
            JTextArea msg = new JTextArea("No reports available");
            scrollPane = new JScrollPane(msg);
        }

        contentPane.add(scrollPane, new GridBagConstraints(0,0,2,1,1,1,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10,15,10,15),0,0));


        //display the chosen data set name in a box
        JLabel labelSelectedReport = new JLabel("Values to Display: ");
        labelSelectedReport.setPreferredSize(new Dimension(200,20));
        labelSelectedReport.setHorizontalAlignment(SwingConstants.LEFT);
        contentPane.add(labelSelectedReport, new GridBagConstraints(0,1,2,1,0,0,
                    GridBagConstraints.SOUTHWEST, GridBagConstraints.NONE, new Insets(5,5,5,5),0,0) );
        textSelectedReport = new JTextField("no selection");
        //textSelectedReport.setBorder(BorderFactory.createLineBorder(Color.BLACK));
        textSelectedReport.setPreferredSize(new Dimension(200,20));
        textSelectedReport.setHorizontalAlignment(SwingConstants.LEFT);
        //textSelectedReport.setBackground(Color.WHITE);
        contentPane.add(textSelectedReport, new GridBagConstraints(0,2,2,1,0,0,
                    GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(5,5,5,5),0,0));

        //provide ok and cancel buttons
        OKButton = new JButton("OK");
        OKButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                OKButton_action();
            }
        });
        contentPane.add(OKButton, new GridBagConstraints(0,3,1,1,0,0,
                    GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(10,60,20,10),1,0));

        cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_action();
            }
        });
        contentPane.add(cancelButton, new GridBagConstraints(1,3,1,1,0,0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10,10,20,10),1,0));
    }

    public Datum getSelectedBook() {
        return selectedBook;
    }

    private DefaultMutableTreeNode getNode(String nodeName) {
        JTree theTree = theHapFrame.getDataTreePanel().getTree();
        TreeModel theModel = theTree.getModel();
        DefaultMutableTreeNode theRoot = (DefaultMutableTreeNode) theModel.getRoot();
        Enumeration enumTree = theRoot.depthFirstEnumeration();
        while (enumTree.hasMoreElements()) {
             DefaultMutableTreeNode node = (DefaultMutableTreeNode) enumTree.nextElement();
             if (!node.isLeaf()) {
                 Datum book = (Datum) node.getUserObject();
                if (book.getName().equals(nodeName)) {
                    return node;
                }
             }
        }
        return null;
    }

    public void OKButton_action() {
        okButtonPressed = true;
        this.setVisible(false);
    }

    public void cancelButton_action() {
        okButtonPressed = false;
        this.setVisible(false);
    }

    public void reportList_ValueChanged() {
        selectedBook = (Datum) reportList.getSelectedValue();
        if (selectedBook.getData() instanceof ReportWriter) {
            ReportWriter aReport = (ReportWriter) selectedBook.getData();
            Object[] columnNames = aReport.getTableColumnNames();
            if(columnNames[0].equals("Haplotype")) {
                textSelectedReport.setText(selectedBook.getName());
                return;
            }
        }
        textSelectedReport.setText("Not a valid haplotype report");
        selectedBook = null;
    }

    public boolean wasOkButtonPressed() {
        return okButtonPressed;
    }


}
