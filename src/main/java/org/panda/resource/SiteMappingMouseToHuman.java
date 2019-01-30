package org.panda.resource;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SiteMappingMouseToHuman extends FileServer
{

	Map<String, Set<String>> mouseToHumanUP;

	Map<String, Map<String, String>> mouseToHumanSiteMap;

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



	public Set<String> mapToHumanProt(String mUP)
	{
		if (mouseToHumanUP.containsKey(mUP)) return mouseToHumanUP.get(mUP);
		return Collections.emptySet();
	}

	public Map<String, List<String>> mapToHumanSite(String mUP, String... sites)
	{
		Map<String, List<String>> result = new HashMap<>();

		if (mouseToHumanUP.containsKey(mUP))
		{
			for (String hUP : mouseToHumanUP.get(mUP))
			{
				String key = mUP + " " + hUP;
				List<String> siteList = new ArrayList<>();
				Map<String, String> siteMap = mouseToHumanSiteMap.get(key);

				for (String mSite : sites)
				{
					String hSite = siteMap.get(mSite);

					if (hSite != null)
					{
						siteList.add(hSite);
					}
				}
				if (!siteList.isEmpty())
				{
					result.put(hUP, siteList);
				}
			}
			return result;
		}
		else return Collections.emptyMap();
	}

	@Override
	public boolean load() throws IOException
	{
		mouseToHumanUP = new HashMap<>();
		mouseToHumanSiteMap = new HashMap<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String hUP = t[2];
			String mUP = t[3];
			String hSite = t[4];
			String mSite = t[5];

			if (mSite.equals("ND")) return;

			double eValue = t[6].equals("--") ? 1 : Double.valueOf(t[6]);
			int mismatch = Integer.valueOf(t[7]);

			if (eValue > 0.01 && mismatch > 3) return;

			if (!mouseToHumanUP.containsKey(mUP)) mouseToHumanUP.put(mUP, new HashSet<>());
			mouseToHumanUP.get(mUP).add(hUP);

			String key = mUP + " " + hUP;

			if (!mouseToHumanSiteMap.containsKey(key)) mouseToHumanSiteMap.put(key, new HashMap<>());
			mouseToHumanSiteMap.get(key).put(mSite, hSite);
		});

		return true;
	}

	public static void main(String[] args)
	{
		Map<String, List<String>> map = get().mapToHumanSite("Q172634", "T23");
		System.out.println("map = " + map);
	}
}
