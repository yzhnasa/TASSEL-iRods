package net.maizegenetics.gwas.modelfitter;

import java.awt.Container;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JRootPane;

import net.maizegenetics.pal.alignment.Alignment;

public class StepwiseOLSModelFitterDialog extends JDialog {
	JCheckBox jboxNested = new JCheckBox("Nest markers within main factor", null, false);
	JList<String> listMainEffects = null;
	int indexOfSelectedEffect = 0;
	
	public StepwiseOLSModelFitterDialog(String[] mainEffects) {
		super();
		
		this.setTitle("Choose Nested Model");
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);
        Container contentPane = getContentPane();
        BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
        contentPane.setLayout(layout);

		int numberOfMainEffects = mainEffects.length;
		if (numberOfMainEffects > 1) {
			
		} else {
			StringBuilder question = new StringBuilder("The dataset contains only one main effect, ");
			question.append(mainEffects[0]);
			contentPane.add(new JLabel(question.toString()));
			contentPane.add(Box.createVerticalStrut(20));
			contentPane.add(jboxNested);
		}
	}
}
