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
        List<Point<T>> Search(double expect, double baseValue);
        List<Point<T>> Search(double expect);
        List<T> SearchContent(double expect, double baseValue);
        List<T> SearchContent(double expect);
        bool Match(double expect, double baseValue);
        bool Match(double expect);
        void Add(Point<T> point);
        double Tolerance();
        ToleranceBy ToleranceType();
        void SetTolerance(double tol);
        void SetToleranceBy(ToleranceBy by);
    }
}
