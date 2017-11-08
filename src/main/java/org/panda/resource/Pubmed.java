package org.panda.resource;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public class Pubmed extends FileServer
{
	Map<String, String> pmcToPm;

	private static Pubmed instance;

	public static Pubmed get()
	{
		if (instance == null) instance = new Pubmed();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"PMC-ids.csv"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/PMC-ids.csv.gz"};
	}

	public String getPMID(String pmcID)
	{
		return pmcToPm.get(pmcID);
	}

	@Override
	public boolean load() throws IOException
	{
		pmcToPm = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0]).skip(1).forEach(line ->
		{
			String[] token = line.split(",");
			if (token.length < 9) return;

			String pmc = token[8];
			String pm = token[9];

			if (pmc.isEmpty() || pm.isEmpty()) return;

			pmcToPm.put(pmc, pm);
		});

		return true;
	}

	public static void main(String[] args)
	{
		String pmc = "PMC3312097";
		System.out.println(pmc + " = " + get().getPMID(pmc));
	}
}
