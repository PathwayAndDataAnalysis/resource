package org.panda.resource;

import org.panda.utility.FileUtil;
import org.panda.utility.SeqAlign;

import java.io.*;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;
import java.util.stream.Collectors;

/**
 * This class is different than other resources. It requires OrthoMam to be downloaded in a local location. It doesn't
 * auto-download it and it doesn't use .panda directory.
 */
public class OrthoMam
{
	private static String ROOT;
	private static Map<String, String> IDMap;

	public static void setRoot(String root)
	{
		ROOT = root;
		IDMap = FileUtil.lines(ROOT + "/ID_map.txt").collect(Collectors.toMap(l -> l.substring(l.indexOf("_")+1, l.indexOf(".")), l -> l.substring(0, l.indexOf("_"))));
	}

	public static Map<String, String> getIDMap()
	{
		return IDMap;
	}

	public static String getAlignmentFilename(String symbol)
	{
		if (!IDMap.containsKey(symbol)) symbol = symbol.toUpperCase();
		if (IDMap.containsKey(symbol))
		{
			return ROOT + "/AA_CDS/" + IDMap.get(symbol) + "_NT_AL_AA.fasta";
		}
		else return null;
	}

	public static String getHumanSeq(String symbol)
	{
		return getSeq(symbol, "Homo_sapiens");
	}

	public static String getSeq(String symbol, String organism)
	{
		return getSeq(symbol, organism, false);
	}

	public static String getSeq(String symbol, String organism, boolean removeGaps)
	{
		String filename = getAlignmentFilename(symbol);
		if (filename == null) return null;
		String seq = "";
		Scanner sc = null;
		try
		{
			sc = new Scanner(new File(filename));
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
			return null;
		}
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith(">") && line.endsWith(organism))
			{
				line = sc.nextLine();
				while (!line.startsWith(">"))
				{
					seq += line;
					if (sc.hasNextLine())
						line = sc.nextLine();
					else break;
				}
			}
		}
		if (removeGaps) seq = seq.replaceAll("-", "");

		if (!seq.isEmpty())
			return seq;
		else return null;
	}



	public static void main(String[] args)
	{
		OrthoMam.setRoot("/home/ozgunbabur/Data/OrthoMaMv12a");
		String seq1 = OrthoMam.getHumanSeq("TP53");
		System.out.println("seq = " + seq1);
		String name = UniProtSequence.get().getNameOfSymbol("TP53", "9606");
		System.out.println("name = " + name);
		String seq2 = UniProtSequence.get().getSequence(name);
		System.out.println("seq = " + seq2);
		System.out.println();
		SeqAlign sa = new SeqAlign(seq1, seq2);
		sa.align();
		sa.print();

		System.out.println();
		sa = new SeqAlign(seq1, "PAPSWPLTTTTSSVPSQ");
//		sa.setFragmentSearch(false);
		sa.align();
		sa.print();

//		int i0 = sa.getMappedIndexS1ToS2(0);
//		int i28 = sa.getMappedIndexS1ToS2(28);
//		int i29 = sa.getMappedIndexS1ToS2(29);
//
//		System.out.println("i0 = " + i0);
//		System.out.println("i28 = " + i28);
//		System.out.println("i29 = " + i29);
	}
}
