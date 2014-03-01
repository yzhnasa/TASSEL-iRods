package net.maizegenetics.analysis.imputation;

import javafx.application.Platform;
import javafx.embed.swing.JFXPanel;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.layout.FlowPane;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;

/**
 * Dialog for driving imputation options for both FILLIN and
 *
 * @author Ed Buckler
 */
public class ImputationDialogPlugin extends AbstractPlugin {
    JFXPanel jfxPanel = new JFXPanel();

    public ImputationDialogPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
        if(isInteractive) {
            JDialog dialog=new JDialog(parentFrame,true);
            dialog.setMinimumSize(new Dimension(400,400));
            initFX(dialog);
            dialog.pack();
            Platform.setImplicitExit(false);
            dialog.setVisible(true);
        }
    }

    private void initFX(JDialog dialog) {
        final JFXPanel fxPanel = new JFXPanel();
        dialog.add(fxPanel);

        Platform.runLater(new Runnable() {
            @Override
            public void run() {
                FlowPane parent = null;
                try {
                    FXMLLoader fxmlLoader = new FXMLLoader(ImputationDialogPlugin.class.getResource("ImputationDialog.fxml"));
                    fxmlLoader.setController(ImputationDialogPlugin.this);
                    parent = (FlowPane)fxmlLoader.load();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                Scene scene = new Scene(parent);
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


    @FXML
    protected void runFILLINHaplotype() {
        System.out.println("The button was clicked!");
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
