package net.maizegenetics.baseplugins.haplotypegraph;

import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.decorators.EdgeWeightLabeller;
import edu.uci.ics.jung.graph.decorators.EdgeWeightLabellerStringer;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import edu.uci.ics.jung.graph.decorators.ToStringLabeller;
import edu.uci.ics.jung.graph.impl.SparseVertex;
import edu.uci.ics.jung.utils.UserData;
import edu.uci.ics.jung.visualization.GraphDraw;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.stats.GLM.ReportWriter;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Aug 12, 2004
 * Time: 4:55:03 PM
 */
public class JungTestPanel extends JPanel {
    //fields

    private final String BLANK = "Blank";
    private final String PHENOTYPE = "Phenotype";
    private final String FREQUENCY = "Frequency";
    private boolean isFrozen;
    private boolean hasPhenotypes = false;
    private HaplotypeFrame theHapFrame;
    JButton freezeButton;
    JButton vertexLabelButton;
    GraphDraw gd;
    Graph g;


    //constructors
    public JungTestPanel() {
        super();
        g = getGraph();
        init();
    }

    public JungTestPanel(Alignment aa, HaplotypeFrame parentFrame) {
        super();
        theHapFrame = parentFrame;
        g = getGraph(aa);
        init();
    }

    //methods
    private void init() {
        //this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        gd = new GraphDraw(g);
        final EdgeWeightLabeller ewl = EdgeWeightLabeller.getLabeller(g);
        EdgeWeightLabellerStringer es = new EdgeWeightLabellerStringer(EdgeWeightLabeller.getLabeller(g));
        edu.uci.ics.jung.visualization.SpringLayout.LengthFunction lengthfunction = new edu.uci.ics.jung.visualization.SpringLayout.LengthFunction() {

            public double getLength(Edge e) {
                return ewl.getWeight(e) * 30;
            }
        };

        //StringLabeller sl = StringLabeller.getLabeller(g);
        StringLabeller sl = ToStringLabeller.setLabellerTo(g, "ToStringLabeller");
        //gd.setGraphLayout(new edu.uci.ics.jung.visualization.SpringLayout(g, f));
        //gd.setGraphLayout(new edu.uci.ics.jung.visualization.FRLayout(g));
        gd.setGraphLayout(new JungSpringLayout(g, lengthfunction));
        gd.setRenderer(new JungHaplotypeRenderer(sl, es));
        gd.addGraphMouseListener(new JungGraphMouseListener((JungHaplotypeGraph) g, gd.getVisualizationViewer()));
        gd.getVisualizationViewer().setBorder(BorderFactory.createLineBorder(Color.BLACK));
        gd.getVisualizationViewer().setPreferredSize(new Dimension(700, 700));
        this.add(gd);

        //add a freeze/unfreeze button to the toolbar
        JPanel toolbarPanel = new JPanel(new GridLayout(5, 1));

        freezeButton = new JButton("Freeze");
        isFrozen = false;
        freezeButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                freezeButton_action();
            }
        });

        //add a button to change the vertex labels
        vertexLabelButton = new JButton(BLANK);
        vertexLabelButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                vertexLabelButton_action();
            }
        });

        //add a button to set the phenotype to be displayed
        JButton setPhenotypeButton = new JButton("Set Phenotype");
        setPhenotypeButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                setPhenotypeButton_action();
            }
        });

        toolbarPanel.add(Box.createVerticalBox());
        toolbarPanel.add(freezeButton);
        toolbarPanel.add(vertexLabelButton);
        //toolbarPanel.add(setPhenotypeButton);
        gd.addTool(toolbarPanel);
    }

    public static Graph getGraph() {
        //PajekNetFile file = new PajekNetFile();
        //Graph g = file.load("\\IdeaProjects\\graph.net");
        //Graph g = TestGraphs.createTestGraph(false);
        File f = new File("C:\\IdeaProjects\\Tassel\\haplotyeTest.txt");
        Graph g = new JungHaplotypeGraph(f);

        return g;
    }

    public static Graph getGraph(Alignment aa) {
        return new JungHaplotypeGraph(aa);
    }

    public void freezeButton_action() {
        if (isFrozen) {
            isFrozen = false;
            freezeButton.setText("Freeze");
            gd.getVisualizationViewer().unsuspend();
        } else {
            isFrozen = true;
            freezeButton.setText("UnFreeze");
            gd.getVisualizationViewer().suspend();
        }
    }

    public void vertexLabelButton_action() {
        //toggle from frequency to phenotype (if exists) to blank
        if (vertexLabelButton.getText().equals(FREQUENCY)) {
            if (hasPhenotypes) {
                vertexLabelButton.setText(PHENOTYPE);
            } else {
                vertexLabelButton.setText(BLANK);
            }
            setAllLabelKeys("freq");

        } else if (vertexLabelButton.getText().equals(PHENOTYPE)) {
            vertexLabelButton.setText(BLANK);
            setAllLabelKeys("phenotype");
        } else {
            vertexLabelButton.setText(FREQUENCY);
            setAllLabelKeys("blank");
        }
        repaint();
    }

    public void setAllLabelKeys(String key) {
        Set vertexSet = g.getVertices();
        Iterator it = vertexSet.iterator();
        while (it.hasNext()) {
            JungHaplotypeVertex jungHaplotypeVertex = (JungHaplotypeVertex) it.next();
            jungHaplotypeVertex.setToStringKey(key);
        }
    }

    public void setPhenotypeButton_action() {
        //get a list of haplotype files and display them for the user to select
        phenotypeChooserDialog pheno = new phenotypeChooserDialog(theHapFrame);
        pheno.show();

        //if the ok button was clicked, change phenotype values to these then display them
        if (pheno.wasOkButtonPressed() && pheno.getSelectedBook() != null) {
            hasPhenotypes = true;
            //create a hashmap of haplotype values keyed on haplotype
            ReportWriter rw = (ReportWriter) pheno.getSelectedBook().getData();
            Object[][] tableData = rw.getTableData();
            HashMap haplotypeValues = new HashMap();
            for (int i = 0; i < tableData.length; i++) {
                haplotypeValues.put(tableData[i][0], tableData[i][1]);
            }

            //add phenotype values to the userDatum for each vertex
            Set vertexSet = g.getVertices();
            Iterator itVertices = vertexSet.iterator();
            while (itVertices.hasNext()) {
                SparseVertex vertex = (SparseVertex) itVertices.next();
                HaplotypeVertex hapVertex = ((JungHaplotypeGraph) g).getHaplotypeVertexFromJungVertex(vertex);
                Object phenotype = haplotypeValues.get(hapVertex.getSequence());
                if (phenotype != null) {
                    vertex.setUserDatum("phenotype", phenotype, UserData.SHARED);
                } else {
                    vertex.setUserDatum("phenotype", "   ", UserData.SHARED);
                }
            }

            //display phenotypes
            setAllLabelKeys("phenotype");
            vertexLabelButton.setText(BLANK);

        }

    }
}

