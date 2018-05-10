package org.panda.resource;

/**
 * @author Ozgun Babur
 */
public class PCPathwayUniprot extends PCPathwayHGNC
{
	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"pc-pathway-uniprot.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://www.pathwaycommons.org/archives/PC2/v9/PathwayCommons9.All.uniprot.gmt.gz"};
	}
}
