package net.maizegenetics.analysis.imputation;

import com.google.common.collect.Range;
import javafx.application.Platform;
import javafx.embed.swing.JFXPanel;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.control.TextField;
import javafx.scene.layout.BorderPane;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import javax.swing.*;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Dialog for driving imputation options for both FILLIN and
 *
 * @author Ed Buckler
 */
public class ImputationDialogPlugin extends AbstractPlugin {
    //JFXPanel jfxPanel = new JFXPanel();
    final JFXPanel fxPanel = new JFXPanel();
    private static Map<TextField, Range> validateFields=new HashMap<>();
    @FXML private static TextField hap_maxOutMiss;
    @FXML private static TextField hap_maxHap;
    @FXML private static TextField hap_hapSize;
    @FXML private static TextField hap_MinSites;
    @FXML private static TextField hap_mxHets;
    @FXML private static TextField hap_mxDiv;

    static {
        validateFields.put(hap_maxOutMiss,Range.closed(0.01,1.0));
        validateFields.put(hap_maxHap,Range.closed(1,10000));
        validateFields.put(hap_hapSize,Range.closed(64,100000000));
        validateFields.put(hap_MinSites,Range.closed(64,100000000));
        validateFields.put(hap_mxHets,Range.closed(0.0,1.0));
        validateFields.put(hap_mxDiv,Range.closed(0.0,1.0));
    }


    JDialog dialog;

    public ImputationDialogPlugin(java.awt.Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
        if(isInteractive) {
            dialog=new JDialog(parentFrame,true);
            dialog.setMinimumSize(new java.awt.Dimension(600,400));
            initFX(dialog);
            dialog.pack();
            Platform.setImplicitExit(false);
            dialog.setVisible(true);
        }
    }

    private void initFX(JDialog dialog) {
        dialog.add(fxPanel);
        final String citationText=this.getCitation();
        Platform.runLater(new Runnable() {
            @Override
            public void run() {
                BorderPane parent = null;
                try {
                    FXMLLoader fxmlLoader = new FXMLLoader(ImputationDialogPlugin.class.getResource("ImputationDialog.fxml"));
                    fxmlLoader.setController(ImputationDialogPlugin.this);
                    parent = (BorderPane)fxmlLoader.load();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                Scene scene = new Scene(parent);
               // citationLabel.setText(citationText);
                fxPanel.setScene(scene);
            }
        });
    }

    public void onAction2(ActionEvent actionEvent) {
        System.out.println("Action2");
    }

    public void hitButton4(ActionEvent actionEvent) {
        System.out.println("hitButton4");

    }


    /**
     * -mxDiv    Maximum divergence from founder haplotype
     + "-mxHet    Maximum heterozygosity of haplotype to even scanned\n"
     + "-hapSize    Preferred haplotype block size in sites (minimum 64); will use the closest multiple of 64 at or below the supplied value\n"
     + "-minPres    Minimum number of present sites within input sequence to do the search\n"
     + "-maxHap    Maximum number of haplotypes per segment\n"
     + "-minTaxa Minimum number of taxa to generate a haplotype\n"
     + "-maxOutMiss  Maximum frequency of missing data in the output haplotype"
     */



    @FXML
    protected void runGUIPlugin() {
        System.out.println("Run button clicked");
        int hapMinSites2=Integer.parseInt(hap_MinSites.getText());
        System.out.println(hapMinSites2);
        dialog.setVisible(false);
    }

    @FXML
    protected void closeGUIPlugin() {
        System.out.println("Close button clicked");
        dialog.setVisible(false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        return null;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return null;
    }

    @Override
    public String getToolTipText() {
        return null;
    }


}
