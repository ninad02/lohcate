package lohcateEnums;

import java.awt.Color;

/** Pastel Color Palette.  Taken from:
 * http://www.tinygorilla.com/Easter_eggs/PallateHex.html
 * 
 * @author Ninad Dewal
 *
 */
public enum ColorPastel {
	
	// GrayScale HEX
	Black(0x000000, 0, 0, 0),	
	White(0xFFFFFF, 255, 255, 255),	
	Gray_10(0xEBEBEB, 235, 235, 235),	
	Gray_15(0xE1E1E1, 225, 225, 225),	
	Gray_20(0xD7D7D7, 215, 215, 215),	
	Gray_25(0xD7D7D7, 204, 204, 204),	
	Gray_30(0xC2C2C2, 194, 194, 194),	
	Gray_35(0xB7B7B7, 183, 183, 183),	
	Gray_40(0xACACAC, 172, 172, 172),	
	Gray_45(0xA0A0A0, 161, 161, 161),	
	Gray_50(0x959595, 149, 149, 149),	
	Gray_55(0x898989, 137, 137, 137),	
	Gray_60(0x7D7D7D, 125, 125, 125),	
	Gray_65(0x707070, 112, 112, 112),	
	Gray_70(0x626262, 99, 99, 99),	
	Gray_75(0x555555, 85, 85, 85),	
	Gray_80(0x464646, 70, 70, 70),	
	Gray_85(0x363636, 54, 54, 54),	
	Gray_90(0x262626, 37, 37, 37),	
	Gray_95(0x111111, 17, 17, 17),	
	
	// Primary Color HEX
	RGB_Red(0xFF0000, 255, 0, 0),	
	RGB_Yelllow(0xFFFF00, 255, 255, 0),	
	RGB_Green(0x00FF00, 0, 255, 0),	
	RGB_Cyan(0x00FFFF, 0, 255, 255),	
	RGB_Blue(0x0000FF, 0, 0, 255),	
	RGB_Magenta(0xFF00FF, 255, 0, 255),	
	CMYK_Red(0xED1C24, 237, 28, 36),	
	CMYK_Yellow(0xFFF200, 255, 242, 0),	
	CMYK_Green(0x00A651, 0, 166, 810),	
	CMYK_Cyan(0x00AEEF, 0, 174, 239),	
	CMYK_Blue(0x2E3192, 46, 49, 146),	
	CMYK_Magenta(0xEC008C, 236, 0, 140),	
	
	// Pastel ColorWheel HEX
	Pastel_Red(0xF7977A, 246, 150, 121),	
	Pastel_Red_Orange(0xF9AD81, 249, 173, 129),	
	Pastel_Yellow_Orange(0xFDC68A, 253, 198, 137),	
	Pastel_Yellow(0xFFF79A, 255, 247, 153),	
	Pastel_Pea_Green(0xC4DF9B, 196, 223, 155),	
	Pastel_Yellow_Green(0xA2D39C, 163, 211, 156),	
	Pastel_Green(0x82CA9D, 130, 202, 156),	
	Pastel_Green_Cyan(0x7BCDC8, 122, 204, 200),	
	Pastel_Cyan(0x6ECFF6, 109, 207, 246),	
	Pastel_Cyan_Blue(0x7EA7D8, 125, 167, 217),	
	Pastel_Blue(0x8493CA, 131, 147, 202),	
	Pastel_Blue_Violet(0x8882BE, 135, 129, 189),	
	Pastel_Violet(0xA187BE, 161, 134, 190),	
	Pastel_Violet_Magenta(0xBC8DBF, 189, 140, 191),	
	Pastel_Magenta(0xF49AC2, 244, 154, 193),	
	Pastel_Magenta_Red(0xF6989D, 245, 152, 157),	
	
	// Light ColorWheel HEX
	Light_Red(0xF26C4F, 242, 108, 79),	
	Light_Red_Orange(0xF68E55, 246, 142, 86),	
	Light_Yellow_Orange(0xFBAF5C, 251, 175, 93),	
	Light_Yellow(0xFFF467, 255, 245, 104),	
	Light_Pea_Green(0xACD372, 172, 211, 115),	
	Light_Yellow_Green(0x7CC576, 124, 197, 118),	
	Light_Green(0x3BB878, 60, 184, 120),	
	Light_Green_Cyan(0x1ABBB4, 28, 187, 180),	
	Light_Cyan(0x00BFF3, 0, 191, 243),	
	Light_Cyan_Blue(0x438CCA, 68, 140, 203),	
	Light_Blue(0x5574B9, 86, 116, 185),	
	Light_Blue_Violet(0x605CA8, 96, 92, 168),	
	Light_Violet(0x855FA8, 133, 96, 168),	
	Light_Violet_Magenta(0xA763A8, 168, 100, 168),	
	Light_Magenta(0xF06EA9, 240, 110, 170),	
	Light_Magenta_Red(0xF26D7D, 242, 109, 125),


	// Pure ColorWheel HEX
	Red(0xED1C24, 242, 109, 125),	
	Red_Orange(0xF26522, 242, 101, 34),	
	Yellow_Orange(0xF7941D, 247, 148, 29),	
	Yellow(0xFFF200, 255, 242, 0),	
	Pea_Green(0x8DC73F, 141, 198, 63),	
	Yellow_Green(0x39B54A, 57, 181, 74),	
	Green(0x00A651, 0, 166, 81),	
	Green_Cyan(0x00A99D, 0, 169, 157),	
	Cyan(0x00AEEF, 0, 174, 239),	
	Cyan_Blue(0x0072BC, 0, 114, 88),	
	Blue(0x0054A6, 0, 84, 166),	
	Blue_Violet(0x2E3192, 49, 49, 146),	
	Violet(0x662D91, 46, 49, 146),	
	Violet_Magenta(0x92278F, 146, 39, 1436),	
	Magenta(0xEC008C, 236, 0, 140),	
	Magenta_Red(0xED145B, 237, 20, 91),
	
	
	// Dark ColorWheel HEX
	Dark_Red(0x9E0B0F, 158, 11, 15),	
	Dark_Red_Orange(0xA0410D, 160, 65, 13),	
	Dark_Yellow_Orange(0xA36209, 136, 98, 10),	
	Dark_Yellow(0xABA000, 171, 160, 0),	
	Dark_Pea_Green(0x598527, 89, 133, 39),	
	Dark_Yellow_Green(0x1A7B30, 25, 123, 48),	
	Dark_Green(0x007236, 0, 114, 54),	
	Dark_Green_Cyan(0x00746B, 0, 116, 107),	
	Dark_Cyan(0x0076A3, 0, 118, 163),	
	Dark_Cyan_Blue(0x004B80, 0, 74, 128),	
	Dark_Blue(0x003471, 0, 52, 113),	
	Dark_Blue_Violet(0x1B1464, 27, 20, 100),	
	Dark_Violet(0x440E62, 68, 14, 98),	
	Dark_Violet_Magenta(0x630460, 99, 4, 96),	
	Dark_Magenta(0x9E005D, 158, 0, 93),	
	Dark_Magenta_Red(0x9E0039, 158, 0, 57),	
	
	// Darker ColorWheel HEX
	Darker_Red(0x790000, 121, 0, 0),	
	Darker_Red_Orange(0x7B2E00, 123, 46, 0),	
	Darker_Yellow_Orange(0x7D4900, 125, 73, 0),	
	Darker_Yellow(0x827B00, 130, 123, 0),	
	Darker_Pea_Green(0x406618, 64, 102, 24),	
	Darker_Yellow_Green(0x005E20, 0, 94, 32),	
	Darker_Green(0x005826, 0, 88, 38),	
	Darker_Green_Cyan(0x005952, 0, 89, 82),	
	Darker_Cyan(0x005B7F, 0, 91, 127),	
	Darker_Cyan_Blue(0x003663, 0, 54, 99),	
	Darker_Blue(0x002157, 0, 33, 87),	
	Darker_Blue_Violet(0x0D004C, 13, 0, 76),	
	Darker_Violet(0x32004B, 50, 0, 75),	
	Darker_Violet_Magenta(0x4B0049, 75, 0, 73),	
	Darker_Magenta(0x7B0046, 123, 0, 70),	
	Darker_Magenta_Red(0x7A0026, 122, 70, 38),	
	
	// Browns ColorWheel HEX
	Pale_Cool_Brown(0xC7B299, 199, 178, 156),	
	Light_Cool_Brown(0x998675, 153, 134, 117),	
	Medium_Cool_Brown(0x736357, 115, 99, 87),	
	Dark_Cool_Brown(0x534741, 83, 71, 65),	
	Darker_Cool_Brown(0x37302D, 55, 48, 45),	
	Pale_Warm_Brown(0xC69C6E, 198, 156, 110),	
	Light_Warm_Brown(0xA67C52, 166, 124, 82),	
	Medium_Warm_Brown(0x8C6239, 140, 98, 57),	
	Dark_Warm_Brown(0x754C24, 117, 76, 36),	
	Darker_Warm_Brown(0x603913, 96, 57, 19),	

	;
	
	

	// ===================================
	// Member variables
	// ===================================
	private int mHexCode;
	private float mValueRed;
	private float mValueGreen;
	private float mValueBlue;
	private Color mColor;
	
	private ColorPastel(int hexCode, float valueRed, float valueGreen, float valueBlue) {
		mHexCode    = hexCode;
		mValueRed   = valueRed;
		mValueGreen = valueGreen;
		mValueBlue  = valueBlue;
		mColor = new Color(mHexCode);
	}
	
	public Color getColor() { return mColor; }
}

