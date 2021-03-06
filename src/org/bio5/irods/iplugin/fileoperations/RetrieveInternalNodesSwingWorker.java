package org.bio5.irods.iplugin.fileoperations;

import java.util.ArrayList;
import java.util.List;

import javax.swing.JOptionPane;
import javax.swing.SwingWorker;
import javax.swing.tree.DefaultMutableTreeNode;

import org.apache.log4j.Logger;
import org.bio5.irods.iplugin.bean.IPlugin;
import org.bio5.irods.iplugin.bean.TapasCoreFunctions;
import org.bio5.irods.iplugin.utilities.IrodsUtilities;
import org.irods.jargon.core.connection.SettableJargonProperties;
import org.irods.jargon.core.pub.io.IRODSFile;
import org.irods.jargon.core.pub.io.IRODSFileFactory;
import org.irods.jargon.core.query.CollectionAndDataObjectListingEntry;

public class RetrieveInternalNodesSwingWorker extends
		SwingWorker<Void, Integer> {

	private Object[] pathForInternalFiles;
	private IPlugin irodsImageJ;
	private static IRODSFile iRodsFile;
	List<CollectionAndDataObjectListingEntry> collectionsUnderGivenAbsolutePath = null;

	public RetrieveInternalNodesSwingWorker(Object[] elements,
			IPlugin irodsImageJ) {
		super();
		this.pathForInternalFiles = elements;
		this.irodsImageJ = irodsImageJ;
	}

	/* Logger instantiation */
	static Logger log = Logger.getLogger(RetrieveInternalNodesSwingWorker.class
			.getName());

	@Override
	public Void doInBackground() throws Exception {

		if (null != irodsImageJ) {
			/* Setting jargon properties */
			SettableJargonProperties jp = new SettableJargonProperties();
			log.info("Default threads : " + jp.getMaxParallelThreads());
			jp.setMaxParallelThreads(10);
			log.info("Threads upgraded to : " + jp.getMaxParallelThreads());
			IRODSFileFactory iRODSFileFactory = irodsImageJ
					.getiRODSFileFactory();
			String finalpath = irodsImageJ.getPathTillHome();
			String loopPath="";
			log.info("ChildPath before adding path seperator" +finalpath);

			String path = null;
			if(irodsImageJ.isHomeDirectoryTheRootNode()){
				path =TapasCoreFunctions.getRootDirectoryPath(irodsImageJ);
				log.info("Path till root: "+path);
			}
			if(!irodsImageJ.isHomeDirectoryTheRootNode()){
				path =TapasCoreFunctions.getHomeDirectoryPath(irodsImageJ);
				log.info("Path till home: "+path);
			}

			/*Problem with finalPath - Pending*/
			for (int i = 0; i < pathForInternalFiles.length; i++) {
				loopPath += IrodsUtilities.getPathSeperator()
						+ pathForInternalFiles[i].toString();
			}
			log.info("loop Path: " +loopPath);
			log.info("Final Path: " +path+loopPath);
			log.info("ChildPath" + finalpath);
			finalpath=path+loopPath;

			iRodsFile = iRODSFileFactory.instanceIRODSFile(finalpath);
			collectionsUnderGivenAbsolutePath = FileOperations
					.retrieveCollectionsUnderGivenPath(iRodsFile, irodsImageJ);

		}
		irodsImageJ.setiRodsFile(iRodsFile);
		done();
		return null;
	}

	@Override
	public void done() {
		List<DefaultMutableTreeNode> childNodesListAfterLazyLoading = parseDirectoryContentsUsingList(
				collectionsUnderGivenAbsolutePath, new DefaultMutableTreeNode());
		if (null != childNodesListAfterLazyLoading) {
			if (childNodesListAfterLazyLoading.size() > 0) {
				irodsImageJ
						.setChildNodesListAfterLazyLoading(childNodesListAfterLazyLoading);
			} else {
				log.error("Child Nodes doesn't exist, empty folder!");
				JOptionPane.showMessageDialog(null, "Empty Folder!", "Error",
						JOptionPane.ERROR_MESSAGE);
				log.info("Removing childNodesListAfterLazyLoading from irodsImageJ");
				irodsImageJ.getChildNodesListAfterLazyLoading().clear();
				log.info("size of irodsImageJ.getChildNodesListAfterLazyLoading() "
						+ irodsImageJ.getChildNodesListAfterLazyLoading()
								.size());
			}
		}
	}

	private List<DefaultMutableTreeNode> parseDirectoryContentsUsingList(
			List<CollectionAndDataObjectListingEntry> listOfCollectionsUnderGivenAbsolutePath,
			DefaultMutableTreeNode node) {
		List<DefaultMutableTreeNode> listOfNodes = new ArrayList<DefaultMutableTreeNode>();
		CollectionAndDataObjectListingEntry fileUnderCollectionAndDataObjectListingEntry = null;
		for (int i = 0; i < listOfCollectionsUnderGivenAbsolutePath.size(); i++) {
			fileUnderCollectionAndDataObjectListingEntry = listOfCollectionsUnderGivenAbsolutePath
					.get(i);

			if (!fileUnderCollectionAndDataObjectListingEntry.isCollection()) {
				// System.out.println("File name" +irodsAccountFile.getName()
				// +":" +irodsAccountFile.getAbsolutePath());
				log.info("File name:"
						+ fileUnderCollectionAndDataObjectListingEntry
								.getNodeLabelDisplayValue()
						+ ":"
						+ fileUnderCollectionAndDataObjectListingEntry
								.getFormattedAbsolutePath());
				DefaultMutableTreeNode child = new DefaultMutableTreeNode(
						fileUnderCollectionAndDataObjectListingEntry
								.getNodeLabelDisplayValue(),
						false);
				listOfNodes.add(child);
			}

			if (fileUnderCollectionAndDataObjectListingEntry.isCollection()) {
				// System.out.println("Direc name" +
				// irodsAccountFile.getName());
				log.info("Direc name:"
						+ fileUnderCollectionAndDataObjectListingEntry
								.getNodeLabelDisplayValue());
				DefaultMutableTreeNode child = new DefaultMutableTreeNode(
						fileUnderCollectionAndDataObjectListingEntry
								.getNodeLabelDisplayValue(),
						true);
				// node.add(child);
				listOfNodes.add(child);
				/*
				 * File[] direcFiles=irodsAccountFile.listFiles(); for(int
				 * i=0;i<direcFiles.length;i++){
				 * //System.out.println("File number" +i +"\n depth:"
				 * +direcFiles.length); log.info("File number:" +i +"\t depth:"
				 * +direcFiles.length);
				 * parseDirectoryContentsUsingList(direcFiles[i], child,
				 * irodsAccount); }
				 */
			}
		}
		return listOfNodes;
	}

}
