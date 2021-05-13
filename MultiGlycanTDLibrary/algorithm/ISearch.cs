using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.algorithm
{
    public enum ToleranceBy
    { Dalton, PPM }
    public interface ISearch<T>
    {
        void Init(List<Point<T>> inputs);
        void Init();
        List<T> Search(double expect, double baseValue);
        List<T> Search(double expect);
        bool Match(double expect, double baseValue);
        bool Match(double expect);
        void Add(Point<T> point);
        void SetTolerance(double tol);
        void SetToleranceBy(ToleranceBy by);
    }
}
