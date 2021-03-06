package org.bio5.irods.iplugin.utilities;

public class Constants {
	public static final String NEW_LINE_STRING = "\n";
	public static final String DEFAULT_STORAGE_RESOURCE = "";
	public static String IMAGEJ_CACHE_FOLDER = System.getProperty("user.home")
			+ IrodsUtilities.getPathSeperator() + ".tapas";
	public static String DEFAULT_PATH_SEPERATOR = "/";
	public static String ZONE_IPLANT = "iplant";
	public static String ZONE_SPX = "spx";
	public static final String HOST_IPLANT = "data.iplantcollaborative.org";
	public static final String HOST_SPX = "spxirods.dyndns.org";
	public static String PORT = "1247";
	public static String HOME_STRING = "home";
	public static String FOLDER_SELECTION_CANCELED = "Folder Selection is canceled";
	public static String JTABBEDPANE_SELECTED_TAB_FILE_INFORMATION = "File Information";
	public static boolean JPROGRESS_SET_STRING_PAINTED = true;
	public static boolean SAVE_PANEL_VISIBILITY = true;
	public static boolean IS_HOME_DIRECTORY_THE_ROOT_NODE = false; /*
																	 * true - if
																	 * you want
																	 * to pull
																	 * everything
																	 * from home
																	 * directory
																	 * (This
																	 * includes
																	 * shared
																	 * files
																	 * too).
																	 * False- if
																	 * you want
																	 * to pull
																	 * collections
																	 * from only
																	 * your
																	 * account
																	 */

	/* Enhanced jargon Properties */
	public static boolean COMPUTE_AND_VERIFY_CHECKSUM_AFTER_TRANSFER_OPTION = false;
	public static int MAX_THREADS = 10;
	public static boolean USE_PARALLEL_TRANSFERS_OPTION = true;

	/* Error Strings */
	public static String ERROR_STRING_CONNECTION_REFUSED = "Connection refused";
	public static String ERROR_STRING_UNKNOWN_HOST = "UnknownHostException";

	/* Property file keys */
	public static final String PROPERTY_FILE_NAME = "tapas.properties";
	public static final String PROPERTY_USER_NAME = "login.username";
	public static final String PROPERTY_ZONE_NAME = "login.zone";
	public static final String PROPERTY_HOST_NAME = "login.host";

	/* OBJstat Details */
	public static final String OBJECT_TYPE_COLLECTION = "COLLECTION";
	public static final String OBJECT_TYPE_DATA_OBJECT = "DATA_OBJECT";

}
