package org.panda.resource;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class SiteMappingMouseToHuman extends FileServer
{
	private static Map<String, String> mapping;

	private static SiteMappingMouseToHuman instance;
	public static SiteMappingMouseToHuman get()
	{
		if (instance == null) instance = new SiteMappingMouseToHuman();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"SiteMappingMouseToHuman.txt"};
	}

	public String map(String uniprotID, String site)
	{
		return mapping.get(uniprotID + "-" + site);
	}

	@Override
	public boolean load() throws IOException
	{
		mapping = new HashMap<>();
		Map<String, Set<String>> synonymMap = new HashMap<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String syms = t[0];
			String mouseID = t[2];
			String hSite = t[3];
			String mSite = t[4];

			double eValue = t[5].equals("--") ? 1 : Double.valueOf(t[5]);
			int mismatch = Integer.valueOf(t[6]);

			if (eValue > 0.01 && mismatch > 3) return;

			String otherNames = t[8];

			String[] sym = syms.split("; ");
			String[] syn = otherNames.equals("NA") ? null : otherNames.split("; ", -1);

			assert otherNames.equals("NA") || syn.length == sym.length :
				"syms=" + syms + "\tsyns=" + otherNames;

			for (int i = 0; i < sym.length; i++)
			{
				mapping.put(mouseID + "-" + mSite, sym[i] + "-" + hSite);

				if (syn != null)
				{
					for (String name : syn[i].split(" "))
					{
						if (!name.isEmpty())
						{
							if (!synonymMap.containsKey(sym[i])) synonymMap.put(sym[i], new HashSet<>());
							synonymMap.get(sym[i]).add(syn[i]);
						}
					}
				}
			}
		});

		return true;
	}
}
