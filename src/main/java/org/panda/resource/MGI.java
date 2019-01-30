package org.panda.resource;

import org.panda.utility.statistics.Histogram;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.stream.Stream;

/**
 * This class provides a mapping between HGNC IDs and approved gene symbols.
 *
 * @author Ozgun Babur
 */
public class MGI extends FileServer
{
	private static MGI instance;

	private Map<String, String> sym2up;
	private Map<String, String> up2sym;

	private Map<String, Set<String>> mouse2human;
	private Map<String, Set<String>> human2mouse;

	public static synchronized MGI get()
	{
		if (instance == null) instance = new MGI();
		return instance;
	}

	public String getSymbol(String uniprotID)
	{
		return up2sym.get(uniprotID);
	}

	public String getUniProtID(String symbol)
	{
		return sym2up.get(symbol);
	}

	public Set<String> getCorrespondingHumanSymbols(String sym)
	{
		return mouse2human.getOrDefault(sym, Collections.emptySet());
	}

//	public Set<String> getCorrespondingHumanSymbols(String upOrSym)
//	{
//		String mUP = sym2up.getOrDefault(upOrSym, upOrSym);
//
//		if (mUP != null)
//		{
//			Set<String> ups = SiteMappingMouseToHuman.get().mapToHumanProt(mUP);
//
//			if (ups != null)
//			{
//				Set<String> hSyms = new HashSet<>();
//				for (String hUP : ups)
//				{
//					String sym = HGNC.get().getSymbol(hUP);
//					if (sym != null) hSyms.add(sym);
//				}
//				return hSyms;
//			}
//		}
//
//		return Collections.emptySet();
//	}

	public Set<String> getAllSymbols()
	{
		return sym2up.keySet();
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"mgi.txt", "mgi-human-mapping.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://www.informatics.jax.org/downloads/reports/MRK_SwissProt.rpt",
			"http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"};
	}

	@Override
	public boolean load() throws IOException
	{
		if (load(getResourceAsStream(getLocalFilenames()[0])))
		{
			human2mouse = new HashMap<>();
			mouse2human = new HashMap<>();

			Scanner scanner = new Scanner(new File(locateInBase(getLocalFilenames()[1])));
			int group = 0;
			Set<String> hSet = new HashSet<>();
			Set<String> mSet = new HashSet<>();

			scanner.nextLine();
			while (scanner.hasNextLine())
			{
				String[] t = scanner.nextLine().split("\t");

				int g = Integer.valueOf(t[0]);
				if (g != group)
				{
					if (!hSet.isEmpty() && !mSet.isEmpty())
					{
						for (String h : hSet)
						{
							for (String m : mSet)
							{
								if (!human2mouse.containsKey(h)) human2mouse.put(h, new HashSet<>());
								if (!mouse2human.containsKey(m)) mouse2human.put(m, new HashSet<>());
								human2mouse.get(h).add(m);
								mouse2human.get(m).add(h);
							}
						}
					}

					hSet.clear();
					mSet.clear();
					group = g;
				}

				if (t[1].equals("human")) hSet.add(t[3]);
				else mSet.add(t[3]);
			}
			// no need to handle the last row since we know it has no human corresponding gene
			return true;
		}
		return false;
	}

	public boolean load(Stream<String> resourceStream) throws IOException
	{
		sym2up = new HashMap<>();
		up2sym = new HashMap<>();

		resourceStream.forEach(line -> {
			String[] token = line.split("\t");
			String sym = token[1];
			String up = token[6];
			sym2up.put(sym, up);
			up2sym.put(up, sym);
		});
		return true;
	}


	public static void main(String[] args)
	{
		System.out.println(get().getCorrespondingHumanSymbols("Akt1"));
	}
}
