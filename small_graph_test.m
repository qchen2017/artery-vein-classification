
Gin.node = [];
Gin.link = [];

% 1
iterator = 1;
Gin.node(iterator).idx = [];
Gin.node(iterator).links = [1 2];
Gin.node(iterator).conn = [2 3];
Gin.node(iterator).numLinks = 2;
Gin.node(iterator).comx = 0;
Gin.node(iterator).comy = 0;

% 2
iterator = iterator + 1;
Gin.node(iterator).idx = [];
Gin.node(iterator).links = [1 3 4];
Gin.node(iterator).conn = [1 -1 -1];
Gin.node(iterator).numLinks = 3;
Gin.node(iterator).comx = 0;
Gin.node(iterator).comy = 0;

% 3
iterator = iterator + 1;
Gin.node(iterator).idx = [];
Gin.node(iterator).links = [2 5 6];
Gin.node(iterator).conn = [1 4 5];
Gin.node(iterator).numLinks = 3;
Gin.node(iterator).comx = 0;
Gin.node(iterator).comy = 0;

% 4
iterator = iterator + 1;
Gin.node(iterator).idx = [];
Gin.node(iterator).links = [5 9 10];
Gin.node(iterator).conn = [3 -1 8];
Gin.node(iterator).numLinks = 3;
Gin.node(iterator).comx = 0;
Gin.node(iterator).comy = 0;

% 5
iterator = iterator + 1;
Gin.node(iterator).idx = [];
Gin.node(iterator).links = [6 7 11 12];
Gin.node(iterator).conn = [3 6 8 9];
Gin.node(iterator).numLinks = 4;
Gin.node(iterator).comx = 0;
Gin.node(iterator).comy = 0;

% 6
iterator = iterator + 1;
Gin.node(iterator).idx = [];
Gin.node(iterator).links = [7 8];
Gin.node(iterator).conn = [5 7];
Gin.node(iterator).numLinks = 2;
Gin.node(iterator).comx = 0;
Gin.node(iterator).comy = 0;

% 7
iterator = iterator + 1;
Gin.node(iterator).idx = [];
Gin.node(iterator).links = [8 17 18];
Gin.node(iterator).conn = [6 -1 -1];
Gin.node(iterator).numLinks = 3;
Gin.node(iterator).comx = 0;
Gin.node(iterator).comy = 0;

% 8
iterator = iterator + 1;
Gin.node(iterator).idx = [];
Gin.node(iterator).links = [10 11 13 14];
Gin.node(iterator).conn = [4 5 -1 -1];
Gin.node(iterator).numLinks = 4;
Gin.node(iterator).comx = 0;
Gin.node(iterator).comy = 0;

% 9
iterator = iterator + 1;
Gin.node(iterator).idx = [];
Gin.node(iterator).links = [12 15 16];
Gin.node(iterator).conn = [5 -1 -1];
Gin.node(iterator).numLinks = 3;
Gin.node(iterator).comx = 0;
Gin.node(iterator).comy = 0;

% --------------------------------------------------------

% 1
iterator = 1;
Gin.link(iterator).n1 = 1;
Gin.link(iterator).n2 = 2;
Gin.link(iterator).point = [];

% 2
iterator = iterator + 1;
Gin.link(iterator).n1 = 1;
Gin.link(iterator).n2 = 3;
Gin.link(iterator).point = [];

% 3
iterator = iterator + 1;
Gin.link(iterator).n1 = 2;
Gin.link(iterator).n2 = -1;
Gin.link(iterator).point = [];

% 4
iterator = iterator + 1;
Gin.link(iterator).n1 = 2;
Gin.link(iterator).n2 = -1;
Gin.link(iterator).point = [];

% 5
iterator = iterator + 1;
Gin.link(iterator).n1 = 3;
Gin.link(iterator).n2 = 4;
Gin.link(iterator).point = [];

% 6
iterator = iterator + 1;
Gin.link(iterator).n1 = 3;
Gin.link(iterator).n2 = 5;
Gin.link(iterator).point = [];

% 7
iterator = iterator + 1;
Gin.link(iterator).n1 = 6;
Gin.link(iterator).n2 = 5;
Gin.link(iterator).point = [];

% 8
iterator = iterator + 1;
Gin.link(iterator).n1 = 6;
Gin.link(iterator).n2 = 7;
Gin.link(iterator).point = [];

% 9
iterator = iterator + 1;
Gin.link(iterator).n1 = 4;
Gin.link(iterator).n2 = -1;
Gin.link(iterator).point = [];

% 10
iterator = iterator + 1;
Gin.link(iterator).n1 = 4;
Gin.link(iterator).n2 = 8;
Gin.link(iterator).point = [];

% 11
iterator = iterator + 1;
Gin.link(1).n1 = 5;
Gin.link(1).n2 = 8;
Gin.link(1).point = [];

% 12
iterator = iterator + 1;
Gin.link(iterator).n1 = 5;
Gin.link(iterator).n2 = 9;
Gin.link(iterator).point = [];

% 13
iterator = iterator + 1;
Gin.link(iterator).n1 = 8;
Gin.link(iterator).n2 = -1;
Gin.link(iterator).point = [];

% 14
iterator = iterator + 1;
Gin.link(iterator).n1 = 8;
Gin.link(iterator).n2 = -1;
Gin.link(iterator).point = [];

% 15
iterator = iterator + 1;
Gin.link(iterator).n1 = 9;
Gin.link(iterator).n2 = -1;
Gin.link(iterator).point = [];

% 16
iterator = iterator + 1;
Gin.link(iterator).n1 = 9;
Gin.link(iterator).n2 = -1;
Gin.link(iterator).point = [];

% 17
iterator = iterator + 1;
Gin.link(iterator).n1 = 7;
Gin.link(iterator).n2 = -1;
Gin.link(iterator).point = [];

% 18
iterator = iterator + 1;
Gin.link(iterator).n1 = 7;
Gin.link(iterator).n2 = -1;
Gin.link(iterator).point = [];



[Gout] = generateGraphOfSegments(Gin);

