"""
answer_writer.py

"""

from io import TextIOWrapper
from typing import List, Tuple, Union


class AnswerWriter:
    """Automates the recording of answers in text files"""

    def __init__(
        self,
        obj: object,
        num_questions: int,
        path: str,
    ) -> None:
        self.obj = obj
        self.type_and_chapter = obj.type_and_chapter
        self.num_questions = num_questions
        self.path = path

    def _multi_part_answer(
        self, answers: List[Tuple[Union[float, str], str, int]]
    ) -> str:
        """Format the answers of multi-part problems"""
        arr = []
        for answer in answers:
            ans, units, sig_figs = answer[0], answer[1], answer[2]
            if isinstance(answer, str):
                arr.append(answer)
            else:
                arr.append(f"{ans:.{sig_figs}} {units}")
        return ", ".join(arr)

    def _format_and_write(
        self,
        question: int,
        answer: Tuple[Union[float, str], str, int],
        file: TextIOWrapper,
    ) -> None:
        """Format single-part answer"""
        q_type, ch_num = self.type_and_chapter[0], self.type_and_chapter[1]
        ans, units, sig_figs = answer[0], answer[1], answer[2]
        if isinstance(answer, str):
            file.write(f"{q_type} {ch_num}.{question}: {answer}\n")
        elif units == "":
            file.write(f"{q_type} {ch_num}.{question}: {ans:.{sig_figs}g}\n")
        else:
            file.write(f"{q_type} {ch_num}.{question}: {ans:.{sig_figs}g} {units}\n")

    def _generate_answers(self) -> List[Tuple[Union[float, str], str, int]]:
        """Iteratively calls solution methods and returns a list of answers"""
        q_type, ch_num = self.type_and_chapter[0], self.type_and_chapter[1]
        n = self.num_questions
        method_names = [f"{q_type.lower()}_{ch_num}_{i}" for i in range(1, n + 1)]
        return [getattr(self.obj, m)() for m in method_names]

    def write_answers(
        self,
    ) -> None:
        """Writes answers to text files"""
        q_type, ch_num = self.type_and_chapter[0], self.type_and_chapter[1]
        answers = self._generate_answers()
        with open(self.path, "a", encoding="utf-8") as file:
            for j, answer in enumerate(answers, start=1):
                if isinstance(answer, List):
                    formatted_ans = self._multi_part_answer(answer)
                    file.write(f"{q_type} {ch_num}.{j}: {formatted_ans}\n")
                else:
                    self._format_and_write(j, answer, file)
            file.write("\n")
