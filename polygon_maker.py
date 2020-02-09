import pyproj
import math
import pprint

class LatLon:
    '''
    緯度経度（世界測地系,度分秒表記）を表すクラス。
    '''

    # 緯度1度あたりの距離[m]。
    DISTANCE_PER_LATITUDE = 111000

    # 関東における経度1度あたりの距離[m]。
    DISTANCE_PER_LONGITUDE = 90163.29	

    def __init__(self, lat, lon):
        '''
        Parameters
        ----------
        lat : int
        　経度。
        lon : int
        　緯度。
        '''
        self.lat = lat
        self.lon = lon
    
    def slide(self, slide_lat, slide_lon):
        '''
        自身を平行移動したLonLatオブジェクトを返す。

        Parameter
        ---------
        slide_lat : int
        　経度方向の平行移動距離[m]。(南へ平行移動する場合は負値)
        slide_lon : int
        　緯度方向の平行移動距離[m]。(西へ平行移動する場合は負値)

        Returns
        -------
        slided_latlon : LatLon
        　平行移動後のLonLatオブジェクト。
        '''

        lat = self.lat + (slide_lat / LatLon.DISTANCE_PER_LATITUDE)
        lon = self.lon + (slide_lon / LatLon.DISTANCE_PER_LONGITUDE) 
        
        return LatLon(lat, lon)
    
    @classmethod
    def calc_distance(cls, a, b):
        '''
        2つの緯度経度オブジェクトの距離[m]を返す。

        Parameters
        ----------
        a,b : LatLon
        　対象となる緯度経度オブジェクト。

        Returns
        -------
        distance : int
        　aとbの距離[m]。
        '''
        grs80 = pyproj.Geod(ellps='GRS80')
        dummy1, dummy2, distance = grs80.inv(a.lon, a.lat, b.lon, b.lat)
        return distance

    @classmethod
    def coord_to_latlon(cls, coord):
        '''
        Coordオブジェクトを緯度経度（世界測地系,度分秒表記）の点に変換する。

        Parameter
        ---------
        coord : Coord
        　変換対象のCoordオブジェクト
        
        Returns
        -------
        　Coordオブジェクトを緯度経度（世界測地系,度分秒表記）に変換したLonLatオブジェクト。
        '''
        return Coord.mid_point.slide(coord.y, coord.x)

    @staticmethod
    def calc_mid_point(a, b):
        '''
        2つの緯度経度オブジェクトの中点を返す。

        Parameters
        ----------
        a,b : LatLon
        　対象となる緯度経度オブジェクト

        Returns
        　aとbの中間地点のLonLatオブジェクトを返す。
        -------
        '''
        return LatLon((a.lat + b.lat)/2, (a.lon + b.lon)/2)

    
class Coord:
    '''
    X-Y座標系での座標を表すクラス。

    Parameters
    ----------
    focus_point_length : int
    　原点から焦点までの距離。
    sin : int
    　正弦値[rad]。
    cos : int
    　余弦値[rad]。
    mid_point : LatLon
    　出発地と目的地の中点を表す LatLon オブジェクト
    '''
    focus_point_length = 0
    sin = 0.0
    cos = 0.0
    mid_point = None

    def __init__(self, x, y):
        '''
        Parameters
        ----------
        x : int
        　x座標[m]。
        y : int
        　y座標[m]。
        '''
        self.x = x
        self.y = y
    
    @classmethod
    def set_class_variables(cls, org, dst):
        '''
        クラス変数をセットする。

        Parameters
        ----------
        '''
        # 焦点距離
        cls.focus_point_length = LatLon.calc_distance(org,dst)
        
        # 正弦、余弦。(この式は dst が org の南西にある場合のみ有効であることに注意。)
        triangle_point = LatLon(dst.lat, org.lon)
        cls.sin = LatLon.calc_distance(org, triangle_point) / LatLon.calc_distance(org, dst)
        cls.cos = LatLon.calc_distance(dst, triangle_point) / LatLon.calc_distance(org, dst)

        # 中点
        cls.mid_point = LatLon.calc_mid_point(org, dst)
    
    def rotate(self):
        '''
        自身を回転させたCoordオブジェクトを返す。
        
        Returns
        -------
        rotated_coord : Coord
        　sin,cosから生成された回転行列にしたがって自信を回転させたCoordオブジェクト。
        '''
        rotated_x = self.x * Coord.cos - self.y * Coord.sin
        rotated_y = self.x * Coord.sin + self.y + Coord.cos
        return Coord(rotated_x, rotated_y)


class EllipicPolygonMaker:
    # 1分間に移動可能な距離[m]
    MOVEBLE_DISTANCE_PER_MIN = 100

    @staticmethod
    def make_ellipic_polygon(org, dst, total_time):
        '''
        出発地・目的地の緯度経度を受け取り、ポリゴンの座標リストを返す。

        Parameters
        ----------
        org : LatLon
        　出発地の緯度経度。
        dst : LatLon
        　目的地の緯度経度。
        total_time : int
        　出発地 -> 経由地 -> 目的地の移動で費やして良い合計時間[min]。

        Returns
        -------
        ellipic_polygon_list : array-like(LatLon)
        　ポリゴンをあらわすLonLatのリスト。

        Notes
        -----
    　　すべての緯度経度は日本測地系度分秒表記。
    　　出発地が目的地よりも北東にある場合のみ動作を保証する。
        '''
        # Coordクラスを初期化
        Coord.set_class_variables(org, dst)
        
        # 移動可能な時間から距離に変換
        detour_distance = EllipicPolygonMaker.MOVEBLE_DISTANCE_PER_MIN * total_time

        # 回転前のポリゴンを取得
        polygon_without_rotate = EllipicPolygonMaker._polygon_without_rotation(Coord.focus_point_length, detour_distance)

        # ポリゴンを回転させ、Coordの世界からLatLonの世界に変換する
        ellipic_polygon_list = []
        for coord in polygon_without_rotate:
            rotated_coord = coord.rotate()
            rotated_latlon = LatLon.coord_to_latlon(rotated_coord)
            ellipic_polygon_list.append(rotated_latlon)
        
        return ellipic_polygon_list

    @staticmethod
    def _obal_point_y(focus_point_length, detour_destance, x):
        '''
        楕円とx座標が与えられた際に、y座標を計算する。
        '''
        return math.sqrt(((detour_destance/2)**2 - focus_point_length**2)*(1 - x**2/(detour_destance/2)**2))

    @staticmethod
    def _obal_point_x(focus_point_length, detour_destance, y):
        '''
        楕円とy座標が与えられた際に、x座標を計算する。
        '''
        return math.sqrt((detour_destance/2)**2 * (1 - y**2/((detour_destance/2)**2 - focus_point_length**2)))

    @staticmethod
    def _polygon_without_rotation(focus_point_length, detour_destance):
        '''
        回転前の座標(x,y)を取得する。

        Parameters:
        -----------
        focus_point_length : int
        　焦点距離[m]。
        detour_distance : int
        　迂回距離[m]。

        Returns:
        --------
        ret_list : array-like (Coord)
        　回転前のポリゴンを表すCoordオブジェクトのリスト。
        '''

        num_points = 8
        px = [-1] * num_points
        py = [-1] * num_points
        # 第一象限(0,1,2)
        px[0] = EllipicPolygonMaker._obal_point_x(focus_point_length, detour_destance, 0)
        py[0] = 0
        px[1] = EllipicPolygonMaker._obal_point_x(focus_point_length, detour_destance, 0)/2
        py[1] = EllipicPolygonMaker._obal_point_y(focus_point_length, detour_destance, px[1])
        px[2] = 0
        py[2] = EllipicPolygonMaker._obal_point_y(focus_point_length, detour_destance, 0)
        # 第二象限(3,4)
        px[3] = px[1] * (-1)
        py[3] = py[1]
        px[4] = px[0] * (-1)
        py[4] = py[0]
        # 第三象限(5,6)
        px[5] = px[3]
        py[5] = py[3] * (-1)
        px[6] = px[2]
        py[6] = py[2] * (-1)
        # 第四象限(7)
        px[7] = px[1]
        py[7] = py[1] * (-1)
        ret_list = []
        for i in range(num_points):
            ret_list.append(Coord(px[i], py[i]))
        return ret_list


if __name__ == "__main__":
    org = LatLon(35.690988, 139.820423)
    dst = LatLon(35.685739, 139.812724)
    total_time = 20 # min
    polygon = EllipicPolygonMaker.make_ellipic_polygon(org, dst, total_time)
    for point in polygon:
        print(str(point.lat) + "," + str(point.lon))
