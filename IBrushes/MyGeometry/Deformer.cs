using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
using System.Text;

namespace MyGeometry
{
	[TypeConverterAttribute(typeof(DeformerConverter))]
	public interface Deformer
	{
		void Deform();
		void Update();
		void Display();
		void Move();
		void MouseDown();
		void MouseUp();
	}

	public class DeformerConverter : ExpandableObjectConverter
	{
		public override bool CanConvertTo(ITypeDescriptorContext context,
										System.Type destinationType)
		{
			if (destinationType == typeof(Deformer))
				return true;

			return base.CanConvertTo(context, destinationType);
		}

		public override object ConvertTo(ITypeDescriptorContext context,
										CultureInfo culture,
										object value,
										System.Type destinationType)
		{
			if (destinationType == typeof(System.String) &&
				 value is Deformer)
			{

				Deformer deformer = (Deformer)value;

				return deformer.ToString();
			}
			return base.ConvertTo(context, culture, value, destinationType);
		}
	}
}
