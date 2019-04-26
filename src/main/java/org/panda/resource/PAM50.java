package org.panda.resource;

import org.panda.utility.ArrayUtil;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Provides genes in OncoKB.
 *
 * @author Ozgun Babur
 */
public class PAM50 extends FileServer
{
	private static PAM50 instance;

	private Map<String, Map<String, Double>> data;
	private Set<String> genes;

	public static synchronized PAM50 get()
	{
		if (instance == null) instance = new PAM50();
		return instance;
	}

	public Set<String> getAllSymbols()
	{
		return genes;
	}

	public Set<String> getUpInBasal()
	{
		return getUpInType("Basal");
	}

	public Set<String> getUpInType(String subtype)
	{
		Map<String, Double> tMap = data.get(subtype);

		return tMap.keySet().stream().filter(g -> data.keySet().stream().filter(name -> !name.equals("Normal"))
			.noneMatch(type -> data.get(type).get(g) > tMap.get(g))).collect(Collectors.toSet());
	}

	public Set<String> getDownInBasal()
	{
		return getDownInType("Basal");
	}

	public Set<String> getDownInType(String subtype)
	{
		Map<String, Double> bMap = data.get(subtype);

		return bMap.keySet().stream().filter(g -> data.keySet().stream().filter(name -> !name.equals("Normal"))
			.noneMatch(type -> data.get(type).get(g) < bMap.get(g))).collect(Collectors.toSet());
	}


	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"PAM.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"https://genome.unc.edu/pubsup/breastGEO/pam50_centroids.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		data = new HashMap<>();
		genes = new HashSet<>();

		String[] header = getResourceAsStream(getLocalFilenames()[0]).findFirst().get().split("\t");

		getResourceAsStream(getLocalFilenames()[0]).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			for (int i = 1; i < t.length; i++)
			{
				if (!data.containsKey(header[i])) data.put(header[i], new HashMap<>());
				data.get(header[i]).put(t[0], Double.valueOf(t[i]));
			}

			genes.add(t[0]);
		});

		return true;
	}

	public static void main(String[] args) throws IOException
	{
		Set<String> downInBasal = get().getDownInBasal();
		Set<String> upInBasal = get().getUpInBasal();

		System.out.println("upInBasal = " + upInBasal);
		System.out.println("downInBasal = " + downInBasal);
	}
}
