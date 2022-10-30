#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace Geometry {

class Vector;

class IShape;
class Point;
class Segment;
class Ray;
class Line;
class Circle;
class Polygon;

class Vector {
 public:
  // Конструтор по умолчанию
  Vector() {
    r_vect_[0] = 0;
    r_vect_[1] = 1;
  }

  Vector(int64_t x_coord, int64_t y_coord) {
    r_vect_[0] = x_coord;
    r_vect_[1] = y_coord;
  }

  Vector(int64_t x_1, int64_t y_1, int64_t x_2, int64_t y_2)
      : Vector(x_2 - x_1, y_2 - y_1) {}

  Vector(const Point& x_coord, const Point& y_coord);
  Vector(const Vector& other) {
    r_vect_[0] = other.r_vect_[0];
    r_vect_[1] = other.r_vect_[1];
  }

  Vector& operator=(const Vector& other) {
    if (this != &other) {
      r_vect_[0] = other.r_vect_[0];
      r_vect_[1] = other.r_vect_[1];
    }
    return *this;
  }

  int64_t& operator[](const int64_t kIndex) {
    if (kIndex == 1) {
      return r_vect_[1];
    }
    return r_vect_[0];
  }

  const int64_t& operator[](const int64_t kIndex) const {
    if (kIndex == 1) {
      return r_vect_[1];
    }
    return r_vect_[0];
  }

  double VecLen() {
    return sqrt(r_vect_[0] * r_vect_[0] + r_vect_[1] * r_vect_[1]);
  }

  Vector& operator+=(const Vector& other);
  Vector& operator-=(const Vector& other);
  Vector& operator*(int64_t alpha);
  Vector& operator*=(Vector& other);
  Vector& operator-();

 private:
  int64_t r_vect_[2];
};

Vector operator+(const Vector& first, const Vector& second) {
  return Vector(first[0] + second[0], second[1] + first[1]);
}

Vector& Vector::operator-() {
  r_vect_[0] = -r_vect_[0];
  r_vect_[1] = -r_vect_[1];
  return *this;
}

Vector& Vector::operator+=(const Vector& other) {
  *this = *this + other;
  return *this;
}

Vector operator-(const Vector& first, const Vector& second) {
  return Vector(first[0] - second[0], second[1] - first[1]);
}

Vector& Vector::operator-=(const Vector& other) {
  *this = *this - other;
  return *this;
}

int64_t operator*(const Vector& first, const Vector& second) {
  return (first[0] * second[0] + first[1] * second[1]);
}

int64_t operator^(const Vector& first, const Vector& second) {
  return (first[0] * second[1] - first[1] * second[0]);
}

Vector& Vector::operator*(int64_t alpha) {
  r_vect_[0] *= alpha;
  r_vect_[1] *= alpha;
  return *this;
}

class IShape {
 public:
  // void Move(const Vector&) — сдвиг фигуры на соответствующий вектор
  virtual IShape& Move(const Vector& vector) = 0;

  // bool ContainsPoint(const Point&) — проверка, содержит ли фигура точку
  virtual bool ContainsPoint(const Point& point) const = 0;

  // bool CrossSegment(const Segment&) — проверка, пересекает ли отрезок фигуру
  virtual bool CrossesSegment(const Segment& segment) const = 0;

  // IShape* Clone() — вернуть указатель на копию фигуры
  virtual std::unique_ptr<IShape> Clone() const = 0;

  // void ToString() — строковое представление фигуры (формат в примерах)
  virtual std::string ToString() = 0;
  virtual ~IShape() = default;
};

class Point : public IShape {
  int64_t x_coord_;
  int64_t y_coord_;

 public:
  Point() : x_coord_(0), y_coord_(0) {}
  Point(int64_t x_0, int64_t y_0) : x_coord_(x_0), y_coord_(y_0) {}

  Point& operator=(const Point& other) {
    if (this != &other) {
      x_coord_ = other.x_coord_;
      y_coord_ = other.y_coord_;
    }
    return *this;
  }

  bool ContainsPoint(const Point& point) const { return (*this == point); }

  Point& Move(const Vector& vector) {
    x_coord_ += vector[0];
    y_coord_ += vector[1];
    return *this;
  }

  std::string ToString() {
    std::ostringstream os;
    os << "Point(" << x_coord_ << ", " << y_coord_ << ")";
    return os.str();
  }

  std::unique_ptr<IShape> Clone() const {
    auto* point = new Point(x_coord_, y_coord_);
    std::unique_ptr<IShape> ptr(point);
    return ptr;
  }

  int64_t GetX() const { return x_coord_; }
  int64_t GetY() const { return y_coord_; }
  bool CrossesSegment(const Segment& segment) const;
  friend bool operator==(const Point& a, const Point& b);
  friend bool operator<(const Point& a, const Point& b);
  friend bool operator>(const Point& a, const Point& b);
  friend Vector operator-(const Point& first, const Point& second);
};

bool operator==(const Point& a, const Point& b) {
  return ((a.x_coord_ == b.x_coord_) && (a.y_coord_ == b.y_coord_));
}

bool operator<(const Point& a, const Point& b) {
  if (a.x_coord_ == b.x_coord_) {
    return (a.y_coord_ < b.y_coord_);
  }
  return ((a.x_coord_ < b.x_coord_));
}

bool operator>(const Point& a, const Point& b) {
  if (a.x_coord_ == b.x_coord_) {
    return (a.y_coord_ > b.y_coord_);
  }
  return ((a.x_coord_ > b.x_coord_));
}

Vector operator-(const Point& first, const Point& second) {
  return Vector(second, first);
}

Vector::Vector(const Point& x_coord, const Point& y_coord) {
  r_vect_[0] = y_coord.GetX() - x_coord.GetX();
  r_vect_[1] = y_coord.GetY() - x_coord.GetY();
}

class Segment : public IShape {
 public:
  Point start, end;

  Segment(const Point& first, const Point& second)
      : start(first), end(second) {}

  bool ContainsPoint(const Point& point) const {
    Vector point_to_first(start, point);
    Vector point_to_second(end, point);
    return ((point_to_first * point_to_second) <= 0 &&
            (point_to_first ^ point_to_second) == 0);
  }

  IShape& Move(const Vector& vector) {
    start.Move(vector);
    end.Move(vector);
    return *this;
  }

  std::unique_ptr<IShape> Clone() const {
    auto* segment = new Segment(start, end);
    std::unique_ptr<IShape> ptr(segment);
    return ptr;
  }

  std::string ToString() {
    std::ostringstream os;
    os << "Segment(Point(" << start.GetX() << ", " << start.GetY()
       << "), Point(" << end.GetX() << ", " << end.GetY() << "))";
    return os.str();
  }

  bool CrossesSegment(const Segment& segment) const;
};

bool Point::CrossesSegment(const Segment& segment) const {
  Vector segment_ftp(segment.start, *this);
  Vector segment_stp(segment.end, *this);
  return ((segment_ftp ^ segment_stp) == 0 && (segment_ftp * segment_stp) <= 0);
}

class Line : public IShape {
 public:
  Point start, end;
  Vector r_vect;
  std::vector<int64_t> form;

  Line(const Point& point_a, const Point& point_b)
      : start(point_a), end(point_b), r_vect(Vector(point_a, point_b)) {
    form.push_back(r_vect[1]);
    form.push_back(-r_vect[0]);
    form.push_back(-r_vect[1] * start.GetX() + r_vect[0] * start.GetY());
  }

  bool ContainsPoint(const Point& point) const {
    return ((point.GetX() - start.GetX()) * r_vect[1] ==
            (point.GetY() - start.GetY()) * r_vect[0]);
  }

  IShape& Move(const Vector& vector) {
    start.Move(vector);
    end.Move(vector);
    r_vect = Vector(start, end);
    form.clear();
    form.push_back(r_vect[1]);
    form.push_back(-r_vect[0]);
    form.push_back(-r_vect[1] * start.GetX() + r_vect[0] * start.GetY());
    return *this;
  }

  bool CrossesSegment(const Segment& segment) const {
    return (((form[0] * segment.start.GetX() + form[1] * segment.start.GetY() +
              form[2]) <= 0 &&
             (form[0] * segment.end.GetX() + form[1] * segment.end.GetY() +
              form[2]) >= 0) ||
            ((form[0] * segment.start.GetX() + form[1] * segment.start.GetY() +
              form[2]) >= 0 &&
             (form[0] * segment.end.GetX() + form[1] * segment.end.GetY() +
              form[2]) <= 0));
  }

  std::unique_ptr<IShape> Clone() const {
    auto* line = new Line(start, end);
    std::unique_ptr<IShape> ptr(line);
    return ptr;
  }

  std::string ToString() {
    std::ostringstream os;
    os << "Line(" << form[0] << ", " << form[1] << ", " << form[2] << ")";
    return os.str();
  }
};

bool CheckIt(Line line_1, Line line_2) {
  Segment seg_lf(line_1.start, line_2.end);
  Segment seg_ls(line_2.start, line_1.end);
  if (line_2.start > line_1.end) {
    seg_ls.start = line_1.end;
    seg_ls.end = line_2.start;
  }
  if (seg_lf.ContainsPoint(line_1.end) && seg_lf.ContainsPoint(line_2.start)) {
    if (line_2.start.GetX() == line_1.end.GetX()) {
      return (line_2.start.GetY() <= line_1.end.GetY());
    }
    return (line_2.start.GetX() <= line_1.end.GetX());
  }
  if (seg_ls.ContainsPoint(line_2.end) && seg_ls.ContainsPoint(line_1.start)) {
    if (line_1.start.GetX() == line_2.end.GetX()) {
      return (line_1.start.GetY() <= line_2.end.GetY());
    }
    return (line_1.start.GetX() <= line_2.end.GetX());
  }
  return true;
}

bool Segment::CrossesSegment(const Segment& segment) const {
  bool flag;
  Segment first_seg(start, end);
  Segment second_sec(segment.start, segment.end);
  if (first_seg.start > first_seg.end) {
    Point temp(start.GetX(), start.GetY());
    first_seg.start = first_seg.end;
    first_seg.end = temp;
  }
  if (second_sec.start > second_sec.end) {
    Point temp(second_sec.start.GetX(), second_sec.start.GetY());
    second_sec.start = second_sec.end;
    second_sec.end = temp;
  }
  Line line_1(first_seg.start, first_seg.end);
  Line line_2(second_sec.start, second_sec.end);
  if (line_2.CrossesSegment(first_seg) && line_1.CrossesSegment(second_sec)) {
    flag = CheckIt(line_1, line_2);
    return flag;
  }
  return false;
}

class Ray : public IShape {
 public:
  Point start, end;
  Vector r_vect;

  Ray(const Point& first, const Point& second)
      : start(first), end(second), r_vect(Vector(first, second)) {}

  bool ContainsPoint(const Point& point) const {
    Vector first_to_point(start, point);
    return (
        (first_to_point[0] * r_vect[1] - first_to_point[1] * r_vect[0] == 0) &&
        (r_vect[0] * first_to_point[0] + r_vect[1] * first_to_point[1] >= 0));
  }

  IShape& Move(const Vector& vector) {
    start.Move(vector);
    return *this;
  }

  bool CrossesSegment(const Segment& segment) const {
    if (ContainsPoint(segment.start) || ContainsPoint(segment.end) ||
        segment.ContainsPoint(start)) {
      return true;
    }
    Vector point_to_sgmt_first(start, segment.start);
    Vector point_to_sgmt_second(start, segment.end);
    if ((point_to_sgmt_first ^ point_to_sgmt_second) == 0) {
      if ((point_to_sgmt_first ^ r_vect) == 0) {
        return ((point_to_sgmt_first * r_vect) >= 0 ||
                (point_to_sgmt_second * r_vect) >= 0);
      }
    }
    return (((point_to_sgmt_first ^ r_vect) > 0 &&
             (r_vect ^ point_to_sgmt_second) > 0 &&
             (point_to_sgmt_first ^ point_to_sgmt_second) > 0) ||
            ((point_to_sgmt_first ^ r_vect) < 0 &&
             (r_vect ^ point_to_sgmt_second) < 0 &&
             (point_to_sgmt_first ^ point_to_sgmt_second) < 0));
  }

  std::unique_ptr<IShape> Clone() const {
    auto* ray = new Ray(start, end);
    std::unique_ptr<Geometry::IShape> ptr(ray);
    return ptr;
  }

  std::string ToString() {
    std::ostringstream os;
    os << "Ray(Point(" << start.GetX() << ", " << start.GetY() << "), Vector("
       << r_vect[0] << ", " << r_vect[1] << "))";
    std::string str = os.str();
    return str;
  }
};

class Polygon : public IShape {
  std::vector<Geometry::Point> polygon_;

 public:
  Polygon(const std::vector<Geometry::Point>& points) : polygon_(points) {}

  static bool Check(Ray ray, Segment side, int64_t count) {
    if (ray.CrossesSegment(side) &&
        !(ray.ContainsPoint(side.start) && ray.ContainsPoint(side.end))) {
      if (ray.ContainsPoint(side.start)) {
        if (side.start.GetY() > side.end.GetY()) {
          count++;
        }
      } else if (ray.ContainsPoint(side.end)) {
        if (side.end.GetY() > side.start.GetY()) {
          count++;
        }
      } else {
        count++;
      }
    }
    return count % 2 != 0;
  }

  bool ContainsPoint(const Point& point) const {
    Point point_another(point.GetX() + 1, point.GetY());
    Ray ray(point, point_another);
    int64_t count = 0;
    for (size_t i = 0; i < polygon_.size() - 1; i++) {
      Segment side(polygon_[i], polygon_[i + 1]);
      if (side.ContainsPoint(point)) {
        return true;
      }
      if (ray.CrossesSegment(side) &&
          !(ray.ContainsPoint(side.start) && ray.ContainsPoint(side.end))) {
        if (ray.ContainsPoint(side.start)) {
          if (side.start.GetY() > side.end.GetY()) {
            count++;
          }
        } else if (ray.ContainsPoint(side.end)) {
          if (side.end.GetY() > side.start.GetY()) {
            count++;
          }
        } else {
          count++;
        }
      }
    }
    Segment side(polygon_[polygon_.size() - 1], polygon_[0]);
    bool res = Check(ray, side, count);
    return res;
  }

  IShape& Move(const Vector& vector) {
    for (size_t i = 0; i < polygon_.size(); i++) {
      polygon_[i].Move(vector);
    }
    return *this;
  }

  bool CrossesSegment(const Segment& segment) const {
    for (size_t i = 0; i < polygon_.size() - 1; i++) {
      Segment side(polygon_[i], polygon_[i + 1]);
      if (segment.CrossesSegment(side)) {
        return true;
      }
    }
    Segment side(polygon_[polygon_.size() - 1], polygon_[0]);
    return (segment.CrossesSegment(side));
  }

  std::unique_ptr<IShape> Clone() const {
    auto* polygon = new Polygon(polygon_);
    std::unique_ptr<Geometry::IShape> ptr(polygon);
    return ptr;
  }

  std::string ToString() {
    std::ostringstream os;
    os << "Polygon(";
    for (size_t i = 0; i < polygon_.size(); i++) {
      if (i == polygon_.size() - 1) {
        os << "Point(" << polygon_[i].GetX() << ", " << polygon_[i].GetY()
           << "))";
      } else {
        os << "Point(" << polygon_[i].GetX() << ", " << polygon_[i].GetY()
           << "), ";
      }
    }
    return os.str();
  }
};

double ToLine(const Point& centre, const Point& first,
              Vector& first_to_second) {
  double result = static_cast<double>(first_to_second[1] * centre.GetX() -
                                      first_to_second[0] * centre.GetY() +
                                      first.GetY() * first_to_second[0] -
                                      first.GetX() * first_to_second[1]) /
                  std::sqrt(first_to_second[1] * first_to_second[1] +
                            first_to_second[0] * first_to_second[0]);
  if (result < 0) {
    return (-result);
  }
  return result;
}

double DoSegment(const Point& centre, const Point& first, const Point& second,
                 Vector& first_to_second, Vector& second_to_first) {
  Vector ftc(first, centre);
  Vector stc(second, centre);
  if (ftc[0] * first_to_second[0] + ftc[1] * first_to_second[1] < 0) {
    return std::sqrt(ftc[0] * ftc[0] + ftc[1] * ftc[1]);
  }
  if (stc[0] * second_to_first[0] + stc[1] * second_to_first[1] < 0) {
    return std::sqrt(stc[0] * stc[0] + stc[1] * stc[1]);
  }
  return ToLine(centre, first, first_to_second);
}

class Circle : public IShape {
  Point centre_;
  int64_t r_;

 public:
  Circle(const Point& centre, int64_t radius) : centre_(centre), r_(radius) {}

  bool ContainsPoint(const Point& point) const {
    Vector dist(centre_, point);
    return (dist.VecLen() <= static_cast<double>(r_));
  }

  IShape& Move(const Vector& vector) {
    centre_.Move(vector);
    return *this;
  }

  bool CrossesSegment(const Segment& segment) const {
    Vector ctf(centre_, segment.start);
    Vector cts(centre_, segment.end);
    if (ctf.VecLen() < static_cast<double>(r_) &&
        cts.VecLen() < static_cast<double>(r_)) {
      return false;
    }
    if ((ctf.VecLen() < static_cast<double>(r_) &&
         cts.VecLen() >= static_cast<double>(r_)) ||
        (ctf.VecLen() >= static_cast<double>(r_) &&
         cts.VecLen() < static_cast<double>(r_))) {
      return true;
    }
    Vector seg_second_to_first(segment.end, segment.start);
    Vector seg_first_to_second(segment.start, segment.end);
    return (DoSegment(centre_, segment.start, segment.end, seg_first_to_second,
                      seg_second_to_first) <= static_cast<double>(r_));
  }

  std::unique_ptr<IShape> Clone() const {
    auto* circle = new Circle(centre_, r_);
    std::unique_ptr<Geometry::IShape> ptr(circle);
    return ptr;
  }

  std::string ToString() {
    std::ostringstream os;
    os << "Circle(Point(" << centre_.GetX() << ", " << centre_.GetY() << "), "
       << r_ << ")";
    return os.str();
  }
};

}  // namespace Geometry

template <class SmartPtrT>
void Delete(const SmartPtrT& ptr) {}

template <class T>
void Delete(T* ptr) {
  delete ptr;
}

void CheckFunctions(const Geometry::IShape* shape_ptr,
                    const Geometry::Point& point_a,
                    const Geometry::Point& point_b) {
  std::cout << "Given shape "
            << (shape_ptr->ContainsPoint(point_a) ? "contains"
                                                  : "does not contain")
            << " point A\n";

  const auto kSegmentAb = Geometry::Segment(point_a, point_b);
  std::cout << "Given shape "
            << (shape_ptr->CrossesSegment(kSegmentAb) ? "crosses"
                                                      : "does not cross")
            << " segment AB\n";

  const auto kVectorAb = point_b - point_a;
  const auto kClonedShapePtr =
      shape_ptr->Clone();  // may return either raw or smart pointer
  std::cout << kClonedShapePtr->Move(kVectorAb).ToString();
}

int main() {
  std::unique_ptr<Geometry::IShape> shape_ptr;

  std::string command;
  std::cin >> command;

  int64_t x = 0;
  int64_t y = 0;
  int64_t a = 0;
  int64_t b = 0;

  if (command == "point") {
    std::cin >> x >> y;
    shape_ptr = std::make_unique<Geometry::Point>(x, y);
  } else if (command == "segment") {
    std::cin >> x >> y >> a >> b;
    shape_ptr = std::make_unique<Geometry::Segment>(Geometry::Point(x, y),
                                                    Geometry::Point(a, b));
  } else if (command == "ray") {
    std::cin >> x >> y >> a >> b;
    shape_ptr = std::make_unique<Geometry::Ray>(Geometry::Point(x, y),
                                                Geometry::Point(a, b));
  } else if (command == "line") {
    std::cin >> x >> y >> a >> b;
    shape_ptr = std::make_unique<Geometry::Line>(Geometry::Point(x, y),
                                                 Geometry::Point(a, b));
  } else if (command == "polygon") {
    size_t n_points = 0;
    std::cin >> n_points;
    std::vector<Geometry::Point> points;
    points.reserve(n_points);
    for (size_t i = 0; i < n_points; ++i) {
      std::cin >> x >> y;
      points.emplace_back(x, y);
    }
    shape_ptr = std::make_unique<Geometry::Polygon>(std::move(points));
  } else if (command == "circle") {
    std::cin >> x >> y;
    const auto kCenter = Geometry::Point(x, y);
    uint64_t radius = 0;
    std::cin >> radius;
    shape_ptr = std::make_unique<Geometry::Circle>(kCenter, radius);
  } else {
    std::cerr << "Undefined command" << std::endl;
    return 1;
  }

  std::cin >> x >> y;
  const auto kPointA = Geometry::Point(x, y);
  std::cin >> x >> y;
  const auto kPointB = Geometry::Point(x, y);

  CheckFunctions(shape_ptr.get(), kPointA, kPointB);
  return 0;
}
