package org.panda.resource.siteeffect;

public enum Feature
{
	PHOSPHORYLATION,
	ACETYLATION,
	METHYLATION,
	UBIQUITINATION,
	GLOBAL_PROTEIN,
	RNA,
	METABOLITE;

	public static Feature getFeat(String code)
	{
		switch (code)
		{
			case "P": return PHOSPHORYLATION;
			case "A": return ACETYLATION;
			case "M": return METHYLATION;
			case "U": return UBIQUITINATION;
			case "G": return GLOBAL_PROTEIN;
			case "R": return RNA;
			case "C": return METABOLITE;
		}
		return null;
	}
}
