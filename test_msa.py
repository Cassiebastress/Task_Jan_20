import unittest
from simplified_msa import (
    aligned_tuples_to_MSA as simplified_aligned_tuples_to_MSA
)


input_data1 = [
    ["a-bcdef", "aAb----"],
    ["a-bcdef", "-Bbc-ef"],
    ["abcde-f-", "---deEfF"]
]

input_data2 = [
    [
        "Today I w-oke up and went to my friend Manu's house for lunch",
        "Today I stood up---------------------------------------------"
    ],
    [
        "Today I woke up and went to my friend Manu--'s house for lunch",
        "----------------------------my friend Cassie's house----------"
    ]
]


class TestTuplesToMSA(unittest.TestCase):
    def test_aligned_tuples_to_MSA_1(self):
        result1 = simplified_aligned_tuples_to_MSA(input_data1)
        expected1 = [
            "a-bcde-f-",
            "aAb------",
            "-Bbc-e-f-",
            "----deEfF"
        ]
        self.assertEqual(result1, expected1)
    
    def test_aligned_tuples_to_MSA_2(self):
        result1 = simplified_aligned_tuples_to_MSA(input_data2)
        expected1 = [
            "Today I w-oke up and went to my friend Manu--'s house for lunch",
            "Today I stood up-----------------------------------------------",
            "-----------------------------my friend Cassie's house----------"
        ]
        self.assertEqual(result1, expected1)

    
if __name__ == '__main__':
    unittest.main()
