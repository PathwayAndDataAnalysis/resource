package org.panda.resource;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.stream.Stream;

/**
 * COSMIC Cancer Gene Census.
 *
 * How to update this resource:
 *
 *       Host name:     sftp-cancer.sanger.ac.uk
 *       Protocol:      sftp
 *       Port:          22
 *       Username:      Your email address
 *       Password:      Your password
 *
 *  Once logged in you will need to download the file from this location /files/grch38/cosmic/v76/
 *  Find the file with the same name, but that have a .csv extension and convert it to a .tsv.
 *
 * @author Ozgun Babur
 */
public class CancerGeneCensus extends FileServer
{
	private static CancerGeneCensus instance;

	private Map<String, String> sym2chr;

	public static synchronized CancerGeneCensus get()
	{
		if (instance == null) instance = new CancerGeneCensus();
		return instance;
	}

	public Set<String> getAllSymbols()
	{
		return sym2chr.keySet();
	}

	public boolean isCancerGene(String sym)
	{
		return sym2chr.containsKey(sym);
	}

	public Set<String> getSymbolsOfChromosome(String chr)
	{
		Set<String> set = new HashSet<>();
		for (String sym : getAllSymbols())
		{
			if (sym2chr.containsKey(sym))
			{
				String c = sym2chr.get(sym);
				String no = c.split("p")[0].split("q")[0];

				String desNo = chr.split("p")[0].split("q")[0];

				if (no.equals(desNo) && c.contains(chr))
				{
					set.add(sym);
				}
			}
		}
		return set;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"cancer_gene_census.tsv"};
	}

	@Override
	public boolean load() throws IOException
	{
		sym2chr = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0]).skip(1).forEach(line -> {
			String[] token = line.split("\t");
			String sym = token[0];
			String chr = token[4];
			sym2chr.put(sym, chr);
		});

		return true;
	}

	public static void main(String[] args)
	{
//		System.out.println(getAllSymbols().contains("LLP"));

		Set<String> genes = HGNC.get().getSymbolsOfChromosome("1q");
		genes.retainAll(CancerGeneCensus.get().getAllSymbols());
		System.out.println("genes = " + genes);
	}
}
