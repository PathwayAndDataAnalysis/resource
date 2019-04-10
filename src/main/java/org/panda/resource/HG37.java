package org.panda.resource;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Stream;

/**
 * This class provides a mapping between HGNC IDs and approved gene symbols.
 *
 * @author Ozgun Babur
 */
public class HG37 extends FileServer
{
	private static HG37 instance;

	private Map<String, int[]> sym2pos;
	private Map<String, String> sym2chr;

	public static synchronized HG37 get()
	{
		if (instance == null) instance = new HG37();
		return instance;
	}

	public Set<String> getAllSymbols()
	{
		return sym2pos.keySet();
	}

	public int[] getLocation(String gene)
	{
		return sym2pos.get(gene);
	}

	public String getChromosome(String gene)
	{
		return sym2chr.get(gene);
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"hg37.gtf"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"};
	}

	@Override
	public boolean load() throws IOException
	{
		sym2pos = new HashMap<>();
		sym2chr = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0]).map(l -> l.split("\t"))
			.filter(t -> t.length >= 9 && t[2].equals("gene"))
			.forEach(t ->
		{
			String gene = t[8].split("; ")[2].split(" ")[1].replaceAll("\"", "");
			String chr = t[0];
			int start = Integer.valueOf(t[3]);
			int end = Integer.valueOf(t[4]);

			sym2pos.put(gene, new int[]{start, end});
			sym2chr.put(gene, chr);
		});

		return true;
	}

	public static void main(String[] args)
	{
		System.out.println("get().getAllSymbols().size() = " + get().getAllSymbols().size());
		System.out.println(Arrays.toString(get().getLocation("TP53")));
	}
}
