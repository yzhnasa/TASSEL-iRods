package org.bio5.irods.iplugin.fileoperations;

import ij.IJ;
import ij.ImagePlus;
import ij.io.Opener;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.SwingWorker;
import javax.swing.tree.TreePath;

import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;
import org.bio5.irods.iplugin.bean.IPlugin;
import org.bio5.irods.iplugin.bean.TapasCoreFunctions;
import org.bio5.irods.iplugin.utilities.Constants;
import org.bio5.irods.iplugin.utilities.IrodsUtilities;
import org.irods.jargon.core.exception.DataNotFoundException;
import org.irods.jargon.core.exception.JargonException;
import org.irods.jargon.core.exception.OverwriteException;
import org.irods.jargon.core.pub.DataObjectAO;
import org.irods.jargon.core.pub.DataTransferOperations;
import org.irods.jargon.core.pub.io.IRODSFile;
import org.irods.jargon.core.pub.io.IRODSFileFactory;
import org.irods.jargon.core.transfer.TransferControlBlock;

public class GetFileFromIrodsSwingWorker extends SwingWorker<Void, Integer> {

	private IRODSFileFactory iRODSFileFactory;
	private String treePath;
	private DataTransferOperations dataTransferOperationsAO;
	private IPlugin iPlugin;
	private DataObjectAO dataObjectAO;
	private TransferControlBlock transferControlBlock;
	IRODSFile sourceIrodsFilePath = null;
	private int executeTimes = 0;

	/* Logger instantiation */
	static Logger log = Logger.getLogger(GetFileFromIrodsSwingWorker.class
			.getName());

	/* Get files from iRODS Server */
	public GetFileFromIrodsSwingWorker(IRODSFileFactory iRODSFileFactory,
			String treePath, IPlugin irodsImagej, JProgressBar progressbar) {
		this.iRODSFileFactory = iRODSFileFactory;
		this.treePath = treePath;
		this.iPlugin = irodsImagej;

	}

	/*
	 * Using SwingWorker-doInBackGround() function to do processing in
	 * background
	 */
	@Override
	public Void doInBackground() throws JargonException {

		log.info("finalTreePath:" + treePath);

		if (null != iPlugin) {
			transferControlBlock = iPlugin.getTransferControlBlock();
			iPlugin.setTransferControlBlock(transferControlBlock);
			iPlugin.setTransferOptions(transferControlBlock
					.getTransferOptions());
			iPlugin.getTransferOptions().setMaxThreads(10);
			dataTransferOperationsAO = iPlugin.getIrodsFileSystem()
					.getIRODSAccessObjectFactory()
					.getDataTransferOperations(iPlugin.getIrodsAccount());
			/*
			 * Check if user requires all files under home directory - this has
			 * performance degradation.
			 */
			if (iPlugin.isHomeDirectoryTheRootNode()) {
				sourceIrodsFilePath =iRODSFileFactory
						.instanceIRODSFile(TapasCoreFunctions.getRootDirectoryPath(iPlugin)+ treePath); 
				log.info("sourceIrodsFilePath" + sourceIrodsFilePath);
						
						/*iRODSFileFactory
						.instanceIRODSFile(IrodsUtilities.getPathSeperator()
								+ iPlugin.getIrodsAccount().getZone()
								+ treePath);*/
			} else {
				sourceIrodsFilePath = iRODSFileFactory
						.instanceIRODSFile(TapasCoreFunctions.getHomeDirectoryPath(iPlugin)+ treePath);
						/*iRODSFileFactory
						.instanceIRODSFile(IrodsUtilities.getPathSeperator()
								+ iPlugin.getIrodsAccount().getZone()
								+ IrodsUtilities.getPathSeperator()
								+ Constants.HOME_STRING + treePath);*/
				log.info("sourceIrodsFilePath" + sourceIrodsFilePath);

			}

			dataObjectAO = iPlugin.getIrodsFileSystem()
					.getIRODSAccessObjectFactory()
					.getDataObjectAO(iPlugin.getIrodsAccount());

			/* Getting MD5 checksum of the current file from iRODS */
			String md5ChecksumLocalFile = null;
			String md5ChecksumServerFile = null;
			try {
				md5ChecksumServerFile = dataObjectAO
						.computeMD5ChecksumOnDataObject(sourceIrodsFilePath);
			} catch (Exception e) {
				log.info("Error while reading MD5 checksum of md5ChecksumServerFile"
						+ e.getMessage());
				JOptionPane
						.showMessageDialog(
								null,
								"Error while reading MD5 checksum of md5ChecksumServerFile!",
								"Error", JOptionPane.ERROR_MESSAGE);
			}
			File destinationLocalFilePath = new File(
					iPlugin.getImageJCacheFolder());
			log.info("sourceIrodsFilePath before inserting file"
					+ sourceIrodsFilePath);
			log.info("destinationLocalFilePath before inserting file"
					+ destinationLocalFilePath);

			/* Getting MD5 checksum of local file, if exists */
			File localFile = new File(
					destinationLocalFilePath.getAbsolutePath()
							+ IrodsUtilities.getPathSeperator()
							+ sourceIrodsFilePath.getName());
			md5ChecksumLocalFile = IrodsUtilities
					.calculateMD5CheckSum(localFile);
			log.info("MD5checksum of iRODS server file: "
					+ md5ChecksumServerFile);
			log.info("MD5checksum of local file: " + md5ChecksumLocalFile);

			if (null != md5ChecksumLocalFile && null != md5ChecksumServerFile
					&& "" != md5ChecksumLocalFile
					&& "" != md5ChecksumServerFile) {
				log.info("MD5 checksum compared - are they Similar files ?"
						+ md5ChecksumLocalFile.equals(md5ChecksumServerFile));

				if (!md5ChecksumLocalFile.equals(md5ChecksumServerFile))
					JOptionPane
							.showMessageDialog(
									null,
									"Local cache directory have files with same name but MD5 checksum is different!",
									"Information",
									JOptionPane.INFORMATION_MESSAGE);
			}
			try {
				if (null != sourceIrodsFilePath) {
					if (null != iPlugin) {

						log.error("Defaulting ErrorWhileUsingGetOperation value to :"
								+ "False");
						iPlugin.setErrorWhileUsingGetOperation(false);

						log.info("Transfer Options in IntraFileStatusCallBack status: "
								+ transferControlBlock.getTransferOptions());
						dataTransferOperationsAO
								.getOperation(
										sourceIrodsFilePath,
										destinationLocalFilePath,
										iPlugin.getIrodsTransferStatusCallbackListener(),
										transferControlBlock);

						if (!iPlugin.isErrorWhileUsingGetOperation()) {
							log.info("Executing openImageUsingImageJ method");
							//openImageUsingImageJ();
							loadDataFileToTassel(executeTimes);
							executeTimes++;
						} else {
							log.error("Error while transfering files");
							JOptionPane
							.showMessageDialog(
									null,
									"Error while transfering files!",
									"Error", JOptionPane.ERROR_MESSAGE);
						}
					}
				}
			} catch (OverwriteException overwriteException) {
				log.error("File with same name already exist in local directory! "
						+ overwriteException.getMessage());
				JOptionPane
						.showMessageDialog(
								null,
								"File with same name already exist in local directory!",
								"Information", JOptionPane.INFORMATION_MESSAGE);

				/* Getting MD5 checksum of local file, if exists */
				File fileInLocal = new File(
						destinationLocalFilePath.getAbsolutePath()
								+ IrodsUtilities.getPathSeperator()
								+ sourceIrodsFilePath.getName());
				md5ChecksumLocalFile = IrodsUtilities
						.calculateMD5CheckSum(fileInLocal);
				log.info("MD5checksum of local file: " + md5ChecksumLocalFile);

				log.info("MD5 checksum compared - Similar files:"
						+ md5ChecksumLocalFile.equals(md5ChecksumServerFile));

				if (!md5ChecksumLocalFile.equals(md5ChecksumServerFile))
					JOptionPane
							.showMessageDialog(null,
									"File names are same but MD5 checksum is different!");

				overwriteException.printStackTrace();
			} catch (DataNotFoundException dataNotFoundException) {
				JOptionPane.showMessageDialog(null, "dataNotFoundException!",
						"Error", JOptionPane.ERROR_MESSAGE);
				log.info("Error while pulling files!"
						+ dataNotFoundException.getMessage());
			} catch (JargonException jargonException) {
				JOptionPane.showMessageDialog(null,
						"Error while pulling files!", "Error",
						JOptionPane.ERROR_MESSAGE);
				log.info("Error while pulling files!"
						+ jargonException.getMessage());
			}

		}
		return null;
	}
	/******************************* Zhong Yang ************************************/
	// When finished download files, then it will load those files to Tassel data tree.
	private void loadDataFileToTassel(int executeTimes) {
		File lopenFile;
		if(null != iPlugin.getImageJCacheFolder()){
			String filePath = iPlugin.getImageJCacheFolder()
					+IrodsUtilities.getPathSeperator()
					+sourceIrodsFilePath.getName();
			File destinationLocalFilePath = new File(
					iPlugin.getImageJCacheFolder());
			iPlugin.getImageJCacheFolder();
			lopenFile = new File(destinationLocalFilePath.getAbsolutePath()
					+ IrodsUtilities.getPathSeperator()
					+ sourceIrodsFilePath.getName());
			iPlugin.getIRodsFileLoadPlugin().setOpenFiles(lopenFile);
			TasselPrefs.putOpenDir(iPlugin.getImageJCacheFolder());
			String myOpenFile = lopenFile.getPath();
			
			DataSet tds = null;
			try {

				if (iPlugin.getIRodsFileLoadPlugin().getFileType() == TasselFileType.Unknown) {
					if (myOpenFile.endsWith(FileLoadPlugin.FILE_EXT_HAPMAP)
							|| myOpenFile
									.endsWith(FileLoadPlugin.FILE_EXT_HAPMAP_GZ)) {
						iPlugin.getIRodsFileLoadPlugin().getLogger().info("guessAtUnknowns: type: "
								+ TasselFileType.Hapmap);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
						tds = iPlugin.getIRodsFileLoadPlugin().processDatum(myOpenFile,
								TasselFileType.Hapmap);
					} else if ((myOpenFile
							.endsWith(FileLoadPlugin.FILE_EXT_TOPM_H5))
							|| (myOpenFile
									.endsWith(FileLoadPlugin.FILE_EXT_TOPM))) {
						iPlugin.getIRodsFileLoadPlugin().getLogger().info("guessAtUnknowns: type: "
								+ TasselFileType.TOPM);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
						tds = iPlugin.getIRodsFileLoadPlugin().processDatum(myOpenFile, TasselFileType.TOPM);
					} else if (myOpenFile
							.endsWith(FileLoadPlugin.FILE_EXT_PLINK_PED)) {
						iPlugin.getIRodsFileLoadPlugin().getLogger().info("guessAtUnknowns: type: "
								+ TasselFileType.Plink);
						String theMapFile = myOpenFile.replaceFirst(
								FileLoadPlugin.FILE_EXT_PLINK_PED,
								FileLoadPlugin.FILE_EXT_PLINK_MAP);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(theMapFile);
						iPlugin.getIRodsFileLoadPlugin().getPlinkLoadPlugin().loadFile(myOpenFile, theMapFile,
								null);
					} else if (myOpenFile
							.endsWith(FileLoadPlugin.FILE_EXT_PLINK_MAP)) {
						iPlugin.getIRodsFileLoadPlugin().getLogger().info("guessAtUnknowns: type: "
								+ TasselFileType.Plink);
						String thePedFile = myOpenFile.replaceFirst(
								FileLoadPlugin.FILE_EXT_PLINK_MAP,
								FileLoadPlugin.FILE_EXT_PLINK_PED);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(thePedFile);
						iPlugin.getIRodsFileLoadPlugin().getPlinkLoadPlugin().loadFile(thePedFile, myOpenFile,
								null);
					} else if (myOpenFile
							.endsWith(FileLoadPlugin.FILE_EXT_SERIAL_GZ)) {
						iPlugin.getIRodsFileLoadPlugin().getLogger().info("guessAtUnknowns: type: "
								+ TasselFileType.Serial);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
						tds = iPlugin.getIRodsFileLoadPlugin().processDatum(myOpenFile,
								TasselFileType.Serial);
					} else if (myOpenFile
							.endsWith(FileLoadPlugin.FILE_EXT_HDF5)) {
						iPlugin.getIRodsFileLoadPlugin().getLogger().info("guessAtUnknowns: type: "
								+ TasselFileType.HDF5);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
						tds = iPlugin.getIRodsFileLoadPlugin().processDatum(myOpenFile, TasselFileType.HDF5);
					} else if (myOpenFile
							.endsWith(FileLoadPlugin.FILE_EXT_VCF)
							|| myOpenFile
									.endsWith(FileLoadPlugin.FILE_EXT_VCF
											+ ".gz")) {
						iPlugin.getIRodsFileLoadPlugin().getLogger().info("guessAtUnknowns: type: "
								+ TasselFileType.VCF);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
						tds = iPlugin.getIRodsFileLoadPlugin().processDatum(myOpenFile, TasselFileType.VCF);
					} else if (myOpenFile
							.endsWith(FileLoadPlugin.FILE_EXT_FASTA)
							|| myOpenFile
									.endsWith(FileLoadPlugin.FILE_EXT_FASTA
											+ ".gz")) {
						iPlugin.getIRodsFileLoadPlugin().getLogger().info("guessAtUnknowns: type: "
								+ TasselFileType.Fasta);
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
						tds = iPlugin.getIRodsFileLoadPlugin().processDatum(myOpenFile, TasselFileType.Fasta);
					} else {
						iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
						tds = iPlugin.getIRodsFileLoadPlugin().guessAtUnknowns(myOpenFile);
					}
				} else {
					iPlugin.getIRodsFileLoadPlugin().setAlreadyLoaded(myOpenFile);
					tds = iPlugin.getIRodsFileLoadPlugin().processDatum(myOpenFile, iPlugin.getIRodsFileLoadPlugin().getFileType());
				}

			} catch (Exception e) {
				e.printStackTrace();
				StringBuilder builder = new StringBuilder();
				builder.append("Error loading: ");
				builder.append(myOpenFile);
				builder.append("\n");
				builder.append(Utils.shortenStrLineLen(
						ExceptionUtils.getExceptionCauses(e), 50));
				String str = builder.toString();
				if (iPlugin.getIRodsFileLoadPlugin().isInteractive()) {
					DialogUtils.showError(str, iPlugin.getIRodsFileLoadPlugin().getParentFrame());
				} else {
					iPlugin.getIRodsFileLoadPlugin().getLogger().error(str);
				}
			}

			if (tds != null) {
				iPlugin.getIRodsFileLoadPlugin().setResult(tds);
				iPlugin.getIRodsFileLoadPlugin().fireDataSetReturned(new PluginEvent(tds, FileLoadPlugin.class));
			}
		}
	}
/***********************************************************************************************/
	@Override
	public void done() {
		/*
		 * Source code in done() method is shifted to openImageUsingImageJ()
		 * method
		 */
	}

	private void openImageUsingImageJ() {

		/* Opening the selected ImageJ */
		Opener imagejOpener = new Opener();
		if (null != iPlugin.getImageJCacheFolder()) {
			String imageFilePath = iPlugin.getImageJCacheFolder()
					+ IrodsUtilities.getPathSeperator()
					+ sourceIrodsFilePath.getName();
			log.info("Current file opened by user: " + imageFilePath);
			ImagePlus imagePlus = imagejOpener.openImage(imageFilePath);
			// ImagePlus imagePlus = IJ.openImage(imageFilePath);

			if (null != imagePlus) {
				iPlugin.setImagePlus(imagePlus);
				log.info("ImagePlus instance is not null and before calling show() function of ImagePlus class");
				imagePlus.show();
				iPlugin.setImageOpened(true);
				log.info("irodsImagej.isImageOpened is set to true");
			} else {
				log.error("ImagePlus instance in GetFileFromIrodsSwingWorker is null and irodsImagej.isImageOpened is false");
				JOptionPane.showMessageDialog(null,
						"File format is not supported by ImageJ!", "Error",
						JOptionPane.ERROR_MESSAGE);
			}
		} else {
			IJ.showMessage("ImageJ is not able to open requested file!");
			IJ.showStatus("ImageJ is not able to open requested file!");
			log.error("ImagePlus instance is null and opening file Failed.");
			JOptionPane.showMessageDialog(null,
					"ImagePlus instance is null and opening file Failed!",
					"Error", JOptionPane.ERROR_MESSAGE);
		}
	}
}
