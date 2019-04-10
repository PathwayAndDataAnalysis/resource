package org.panda.resource;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SiteMappingRatToHuman  extends SiteMappingMouseToHuman
{
	private static SiteMappingRatToHuman instance;

	public static SiteMappingRatToHuman get()
	{
		if (instance == null) instance = new SiteMappingRatToHuman();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"SiteMappingRatToHuman.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		mouseToHumanUP = new HashMap<>();
		mouseToHumanSiteMap = new HashMap<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String hUP = t[1];
			String mUP = t[2];
			String hSite = t[3];
			String mSite = t[4];

			if (mSite.equals("ND")) return;

			double eValue = t[5].equals("--") ? 1 : Double.valueOf(t[5]);
			int mismatch = Integer.valueOf(t[6]);

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
		Map<String, List<String>> map = get().mapToHumanSite("D3Z9J3", "F40");
		System.out.println("map = " + map);
	}
}
