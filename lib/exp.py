class retype_Exception(Exception):
    def __init__(self, tp):
        err = '重复定义了' + tp
        Exception.__init__(self, err)


class repara_Exception(Exception):
    def __init__(self, tp):
        err = '重复定义了' + tp + '的参数'
        Exception.__init__(self, err)
