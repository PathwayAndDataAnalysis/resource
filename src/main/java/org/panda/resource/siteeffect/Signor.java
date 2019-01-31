package org.panda.resource.siteeffect;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;

/**
 * Serves the Signor database.
 *
 * @author Ozgun Babur
 */
public class Signor extends SiteEffectServer
{
	static Signor instance;

	public static Signor get()
	{
		if (instance == null) instance = new Signor();
		return instance;
	}

	public static synchronized void initSingletonEmpty()
	{
		instance = new Signor()
		{
			@Override
			public synchronized boolean init() throws IOException
			{
				return true;
			}
		};

		instance.typeMap = new HashMap<>();
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"Signor_all_data_08_09_17.csv"};
	}

	@Override
	public boolean load() throws IOException
	{
		typeMap = new HashMap<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).skip(1).map(l -> l.split(";"))
			.filter(t -> t.length > 12)
			.filter(t -> t[9].contains("phospho"))
			.filter(t -> t[12].equals("9606"))
			.filter(t -> !t[10].isEmpty())
			.forEach(t ->
		{
			String gene = t[4];
			String effect = t[8];
			String mechanism = t[9];
			String residue = t[10];

			int edgeSign = mechanism.startsWith("de") ? -1 : 1;
			int effSign = effect.startsWith("up") ? 1 : effect.startsWith("down") ? -1 : 0;
			int sign = edgeSign * effSign;

			if (sign == 0) return;

			if (residue.startsWith("Ser")) residue = residue.replace("Ser", "S");
			else if (residue.startsWith("Thr")) residue = residue.replace("Thr", "T");
			else if (residue.startsWith("Tyr")) residue = residue.replace("Tyr", "Y");
			else if (residue.startsWith("His")) residue = residue.replace("His", "H");
			else throw new RuntimeException("Different residue: " + residue);

			if (!typeMap.containsKey(gene)) typeMap.put(gene, new HashMap<>());
			typeMap.get(gene).put(residue, sign);
		});

		return true;
	}

	public static void main(String[] args)
	{
	}
}
