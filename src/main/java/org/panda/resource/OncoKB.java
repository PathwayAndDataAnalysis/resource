package org.panda.resource;

import org.panda.utility.ArrayUtil;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Provides genes in OncoKB.
 *
 * @author Ozgun Babur
 */
public class OncoKB extends FileServer
{
	private static OncoKB instance;

	private Map<String, Map<String, String[]>> data;

	public static synchronized OncoKB get()
	{
		if (instance == null) instance = new OncoKB();
		return instance;
	}

	public Set<String> getAllSymbols()
	{
		return data.keySet();
	}

	public boolean isCancerGene(String sym)
	{
		return data.containsKey(sym);
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"oncoKB.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://oncokb.org/api/v1/utils/allAnnotatedVariants.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		data = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0], StandardCharsets.ISO_8859_1).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (!data.containsKey(t[0])) data.put(t[0], new HashMap<>());

			data.get(t[0]).put(t[1], ArrayUtil.getTail(t, 2));
		});

		return true;
	}

	public void printMatchReport(Map<String, Set<String>> muts)
	{
		muts.keySet().stream().sorted().forEach(gene ->
		{
			if (data.containsKey(gene))
			{
				for (String mut : muts.get(gene))
				{
					if (data.get(gene).containsKey(mut))
					{
						System.out.println(gene + "\t" + mut + "\t" +  Arrays.toString(data.get(gene).get(mut)));
					}
					else if ((mut.endsWith("*") || mut.contains("fs")) &&
						data.get(gene).containsKey("Truncating Mutations"))
					{
						System.out.println(gene + "\t" + mut + " (trunc)"  + "\t" +  Arrays.toString(data.get(gene).get("Truncating Mutations")));
					}
					else
					{
						System.out.println(gene + "\t" + mut);
					}
				}
			}
		});
	}

	public static void main(String[] args) throws IOException
	{
		Map<String, Set<String>> map = new HashMap<>();

		Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/Patient1-revisit/data/met1-muts.maf"))
			.map(l -> l.split("\t")).filter(t -> t.length > 41).filter(t -> !t[41].isEmpty()).forEach(t ->
		{
			if (!map.containsKey(t[0])) map.put(t[0], new HashSet<>());
			map.get(t[0]).add(t[41].substring(2));
		});

		get().printMatchReport(map);
	}
}
