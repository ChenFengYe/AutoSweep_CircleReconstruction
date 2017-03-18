using System;
using System.Collections;
using System.Collections.Generic;

namespace System.Collections.Generic
{
	public class Pair<T1, T2>
	{
		public T1 Item1 { set; get; }
		public T2 Item2 { set; get; }

		public Pair(T1 item1, T2 item2)
		{
			this.Item1 = item1;
			this.Item2 = item2;
		}

		public override bool Equals(object obj)
		{
			Pair<T1, T2> p = obj as Pair<T1, T2>;
			if (p == null) return false;

			if (p.Item1.Equals(this.Item1) && p.Item2.Equals(this.Item2))
				return true;
			else
				return p.Item1.Equals(this.Item2) && p.Item2.Equals(this.Item1);
		}
	}
}
