## Немного информации про структуру репозитория.

---

- `julia_code/`
Здесь хранится код программы. На данный момент всё не очень хорошо структурировано и отлажено, работа ведётся.
- `julia_code/Examples`
Конкретные примеры.  
`example_1.jl` --- Решение первой краевой задачи для одномерного квазилинейного уравнения теплопроводности из статьи Самарского  
`example_2.jl` --- то же самое, что и `example_1.jl`, но решение на неравномерной сетке
`example_3.jl` --- Решение первой краевой задачи для двумерного линейного уравнения теплопроводности (собственный выдуманный пример для проверки)  
`quazi_2.jl` --- Решение первой краевой задачи для двумерного квазилинейного уравнения теплопроводности из статьи Самарского.
- `julia_code/output`  
Вывод всех программ находится в этой папке.


**Видно, что если использовать _статическую_ неравномерную сетку, то порядок сходимости схемы определяется всё-равно _максимальным шагом_ на сетке** (подробнее напишу в отчёте).
- `julia_code/old.ipynb` --- старый файл
- `report/`
Здесь будет делаться отчёт о НИРС. Итоговый документ --- `main.pdf` (его и стоит просматривать).  
- `remarks.md` --- в этом файле можно оставлять любые замечания/предложения по работе.


## Адаптивные сетки
Начал со статического случая. Хочу пока просто создать для примера неконформную сетку и посмотреть, как будет работать алгоритм с "коррекцией потоков"


## Текст диплома
Потихоньку пишу, что могу. Если честно, сложно что-либо писать про Локально-адаптивные сетки пока, т.к. нет результатов:( 
Очень хочу получить хоть какие-то работающие программы, графики, и на основе них уже отталкиваться, что именно писать.
