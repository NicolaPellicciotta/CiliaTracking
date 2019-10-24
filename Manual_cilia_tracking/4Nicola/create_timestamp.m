function [ template ] = create_timestamp( time, flag_whiteonblack )
%create_timestamp creates a tiny matrix with the time, that has to be burned
%into the frames. 
%   Detailed explanation goes here

if nargin < 2 | isempty(flag_whiteonblack)
    flag_whiteonblack = false;
end

my_font = create_font;

ch = 12; %char height
cw = 8; %char_width

% str = sprintf('%08.2f s',time);
str = sprintf('%07.4f s',time);
template = zeros(ch,length(str)*(cw),'uint8');

for cc = 1:length(str)
    
    ind = strcmp(str(cc),{my_font(:).name});
    if sum(ind) == 1
        template(:,cw*(cc-1)+1:cw*cc) = my_font(ind).IM;
    else %put blank
        template(:,cw*(cc-1)+1:cw*cc) = 255*ones(ch,cw,'uint8');
    end
    
end

if flag_whiteonblack
    template = imcomplement(template);
end

end




%% Apologies for the messiness and the hard-codedness


function my_font = create_font()

my_font(1).name = '0';
my_font(1).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     0     0     0     0     0     1     1
    0     0     1     1     1     0     0     1
    0     0     1     1     0     0     0     1
    0     0     1     0     0     0     0     1
    0     0     1     0     1     0     0     1
    0     0     0     0     1     0     0     1
    0     0     0     1     1     0     0     1
    0     0     1     1     1     0     0     1
    1     0     0     0     0     0     1     1
    1     1     1     1     1     1     1     1
    ];


my_font(2).name = '1';
my_font(2).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     0     1     1     1
    1     1     1     0     0     1     1     1
    1     0     0     0     0     1     1     1
    1     1     1     0     0     1     1     1
    1     1     1     0     0     1     1     1
    1     1     1     0     0     1     1     1
    1     1     1     0     0     1     1     1
    1     1     1     0     0     1     1     1
    1     0     0     0     0     0     0     1
    1     1     1     1     1     1     1     1
    ];

my_font(3).name = '2';
my_font(3).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     0     0     0     0     1     1
    1     0     0     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     1     1     1     1     0     0     1
    1     1     1     1     0     0     1     1
    1     1     1     0     0     1     1     1
    1     1     0     0     1     1     1     1
    1     0     0     1     1     0     0     1
    1     0     0     0     0     0     0     1
    1     1     1     1     1     1     1     1
    ];

my_font(4).name = '3';
my_font(4).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     0     0     0     0     1     1
    1     0     0     1     1     0     0     1
    1     1     1     1     1     0     0     1
    1     1     1     1     1     0     0     1
    1     1     1     0     0     0     1     1
    1     1     1     1     1     0     0     1
    1     1     1     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     1     0     0     0     0     1     1
    1     1     1     1     1     1     1     1
    ];

my_font(5).name = '4';
my_font(5).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     0     0     1     1
    1     1     1     0     0     0     1     1
    1     1     0     0     0     0     1     1
    1     0     0     1     0     0     1     1
    0     0     1     1     0     0     1     1
    0     0     0     0     0     0     0     1
    1     1     1     1     0     0     1     1
    1     1     1     1     0     0     1     1
    1     1     1     0     0     0     0     1
    1     1     1     1     1     1     1     1
    ];

my_font(6).name = '5';
my_font(6).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     0     0     0     0     0     0     1
    1     0     0     1     1     1     1     1
    1     0     0     1     1     1     1     1
    1     0     0     0     0     0     1     1
    1     0     0     1     1     0     0     1
    1     1     1     1     1     0     0     1
    1     1     1     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     1     0     0     0     0     1     1
    1     1     1     1     1     1     1     1
    ];

my_font(7).name = '6';
my_font(7).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     0     0     0     1     1
    1     1     0     0     1     1     1     1
    1     0     0     1     1     1     1     1
    1     0     0     1     1     1     1     1
    1     0     0     0     0     0     1     1
    1     0     0     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     1     0     0     0     0     1     1
    1     1     1     1     1     1     1     1
    ];

my_font(8).name = '7';
my_font(8).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    0     0     0     0     0     0     0     1
    0     0     1     1     1     0     0     1
    0     0     1     1     1     0     0     1
    1     1     1     1     1     0     0     1
    1     1     1     1     0     0     1     1
    1     1     1     0     0     1     1     1
    1     0     0     1     1     1     1     1
    1     0     0     1     1     1     1     1
    1     0     0     1     1     1     1     1
    1     1     1     1     1     1     1     1
    ];

my_font(9).name = '8';
my_font(9).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     0     0     0     0     1     1
    1     0     0     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     1     0     0     0     0     1     1
    1     0     0     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     1     0     0     0     0     1     1
    1     1     1     1     1     1     1     1
    ];

my_font(10).name = '9';
my_font(10).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     0     0     0     0     1     1
    1     0     0     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     0     0     1     1     0     0     1
    1     1     0     0     0     0     1     1
    1     1     1     1     0     0     1     1
    1     1     1     1     0     0     1     1
    1     1     1     0     0     1     1     1
    1     1     0     0     0     1     1     1
    1     1     1     1     1     1     1     1
    ];

my_font(11).name = '.';
my_font(11).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     0     0     0     1     1     1
    1     1     0     0     0     1     1     1
    1     1     1     1     1     1     1     1
    ];

my_font(12).name = ':';
my_font(12).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     0     0     0     1     1     1
    1     1     0     0     0     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     0     0     0     1     1     1
    1     1     0     0     0     1     1     1
    1     1     1     1     1     1     1     1
    ];

my_font(13).name = 's';
my_font(13).IM =255.*[...
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1
    1     1     0     0     0     0     1     1
    1     0     0     1     1     0     0     1
    1     1     0     0     1     1     1     1
    1     1     1     1     0     0     1     1
    1     0     0     1     1     0     0     1
    1     1     0     0     0     0     1     1
    1     1     1     1     1     1     1     1
    ];


end