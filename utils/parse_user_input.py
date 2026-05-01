from typing import Optional


def parse_user_input(
        number_parse: str = None,
        prompt: str = "请输入数字",
        max_index: int = None,
        whitelist: str = None
) -> Optional[list[int] | str]:
    """解析输入的数字对应序号

    Args:
        number_parse: 输入数字
        prompt: 输入提示词
        max_index: 最大索引值
        whitelist: 字母白名单

    Returns:
        解析的序号列表
    """
    # 防止用户手误循环判断
    while True:
        # 获取输入并根据，分割、列表化
        parse = input(prompt) if not number_parse else number_parse
        if not parse:
            print("输入不能为空！")
            continue
        # 白名单直接输出
        if whitelist:
            if parse.strip() in whitelist:
                return parse
        parse_split = [p.strip() for p in parse.split(",")]

        # 结果列表和错误捕获初始化
        result_list = []
        error_catcher = False

        # 实现解读“：”
        for index in range(len(parse_split)):
            # 清洗空格
            single_parse = str(parse_split[index]).strip()

            if ":" not in single_parse:
                if not is_digit(list(single_parse)):
                    error_catcher = True
                    break
                result_list.append(int(single_parse))
            else:
                if not is_digit(single_parse.split(":")):
                    error_catcher = True
                    break
                pre, post = single_parse.split(":")
                pre, post = int(pre), int(post)
                result_list.extend(range(pre, post + 1))

        result_list = sorted(list(set(result_list)))

        # 检测是否超过上限
        if not error_catcher and max_index is not None:
            for num in result_list:
                is_out_of_range = True if num > max_index else False
                if is_out_of_range:
                    print(f"最大序号为{max_index},您的输入过大")
                    error_catcher = True
                    break

        # 确认是否选择正确
        if not error_catcher:
            is_confirm = input(f"现在选择的序号是{result_list},是否确认([y]/n):")
            if is_confirm == "n":
                print("根据您的确认，本次输入失效，请重新输入...")
                continue

        # 不合规输入将重新读取,若传入字符串则直接结束循环
        if error_catcher and number_parse:
            print("不合格的字符")
            break
        elif error_catcher and not number_parse:
            print("不合格的字符")
            continue

        return result_list


def is_digit(ls: list) -> bool:
    """简单判断列表里是否都是数字

    Args:
        ls: 输入的列表

    Returns:
        若都是数值返回True否则False
    """
    for num in ls:
        if isinstance(num, int):
            continue
        if not num.isdigit():
            return False
    return True


def main():
    flag = False
    test_whitelist = "m"
    res = parse_user_input(number_parse="a")
    for i in res:
        flag = True if isinstance(i, int) or i in test_whitelist else False
        if not flag:
            break

    if flag:
        print(res)
        print("Done!")


if __name__ == "__main__":
    main()